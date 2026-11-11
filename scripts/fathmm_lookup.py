import math
import re
import sqlite3
import sys
import types
import zipfile


def load_annotator(cursor, code_archive):
    # The CLI base class is irrelevant to the per-variant annotation method.
    # Use the official module code with the already opened read-only database.
    with zipfile.ZipFile(code_archive) as archive:
        source = archive.read("fathmm.py")
    cravat = types.ModuleType("cravat")
    cravat.BaseAnnotator = object
    cravat.InvalidData = ValueError
    previous = sys.modules.get("cravat")
    sys.modules["cravat"] = cravat
    namespace = {"__name__": "fathmm_241_local"}
    try:
        exec(compile(source, "OpenCRAVAT/fathmm/2.4.1/fathmm.py", "exec"), namespace)
    finally:
        if previous is None:
            del sys.modules["cravat"]
        else:
            sys.modules["cravat"] = previous
    annotator = namespace["CravatAnnotator"]()
    annotator.curs = cursor
    return annotator


def parts(value):
    return str(value).split(";") if value is not None else []


def numeric(value):
    if value in {"", "."}:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("Non-finite FATHMM value in source")
    return value


def score_variant(row, annotator, cursor):
    rank, samples, gene, chrom, pos, ref, alt, _, hgvs_p, _, _, transcript = row
    assert re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)", chrom) and alt in "ACGT" and len(alt) == 1
    result = dict(rank=rank, samples=samples, gene=gene, chrom=chrom, pos=pos, ref=ref, alt=alt,
                  hgvs_p=hgvs_p, workbook_transcript=transcript, batch=(rank - 1) // 4 + 1,
                  module_status="NO_SCORE", error=None, scores=None, transcripts=None,
                  predictions=None, min_score=None, rank_score=None, score_source=None,
                  workbook_transcript_score=None, source_score_raw=None, source_transcripts_raw=None,
                  source_predictions_raw=None, source_rank_score_raw=None, source_status="NO_RECORD")
    source_rows = list(dict.fromkeys(cursor.execute(
        f"SELECT ens_tid, ens_pid, fathmm_score, fathmm_rscore, fathmm_pred FROM {chrom} WHERE pos=? AND alt=?",
        (pos, alt),
    )))
    multiple_rows = len(source_rows) > 1
    output = {}
    if multiple_rows:
        result["module_status"] = "MULTIPLE_ROWS"
        target = transcript.split(".")[0] if transcript else None
        exact = [
            candidate for candidate in source_rows
            if target and target in {t.split(".")[0] for t in parts(candidate[0])}
        ]
        if len(exact) != 1:
            result.update(source_status="AMBIGUOUS", error="Multiple database rows with no unique transcript match")
            return result
        source_rows = exact
        result["error"] = "Multiple database rows; selected the unique match to the workbook transcript"
    else:
        try:
            output = annotator.annotate({"chrom": chrom, "pos": pos, "ref_base": ref, "alt_base": alt})
            if output:
                result["module_status"] = "SCORE_RETURNED"
        except (ValueError, TypeError, AttributeError, sqlite3.Error) as exc:
            result.update(module_status="ERROR", error=f"{type(exc).__name__}: {exc}")
    if not source_rows:
        return result
    tids, _, raw, rank_score, predictions = source_rows[0]
    result.update(source_score_raw=raw, source_transcripts_raw=tids,
                  source_predictions_raw=predictions, source_rank_score_raw=rank_score,
                  source_status="NO_NUMERIC_SCORE")
    all_scores, all_tids, all_preds = parts(raw), parts(tids), parts(predictions)
    valid = []
    for i, score in enumerate(all_scores):
        try:
            if numeric(score) is not None:
                valid.append((i, score))
        except (ValueError, TypeError):
            continue
    if not valid:
        return result
    score_source = "OpenCRAVAT FATHMM 2.4.1"
    if multiple_rows:
        score_source = "FATHMM v2.3 database, unique transcript match"
    elif not output:
        score_source = "FATHMM v2.3 database, numeric values recovered"
    result.update(scores=";".join(s for _, s in valid), min_score=min(float(s) for _, s in valid),
                  score_source=score_source, source_status="NUMERIC_SCORE")
    if output:
        result["rank_score"] = output.get("fathmm_rscore")
    if len(all_tids) == len(all_scores):
        result["transcripts"] = ";".join(all_tids[i] for i, _ in valid)
        if transcript:
            target = transcript.split(".")[0]
            matched = {float(s) for i, s in valid if all_tids[i].split(".")[0] == target}
            if len(matched) == 1:
                result["workbook_transcript_score"] = matched.pop()
    if len(all_preds) == len(all_scores):
        result["predictions"] = ";".join("" if all_preds[i] == "." else all_preds[i] for i, _ in valid)
    return result
