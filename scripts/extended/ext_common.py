"""Shared helpers for the extended analyses (no sample identifiers are stored here)."""
import string

import pandas as pd

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
GT_CODE = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2}
HCM_GENES = ["ACTC1", "ALPK3", "CSRP3", "FHOD3", "JPH2", "KLHL24", "MYBPC3", "MYH7", "MYL2", "MYL3", "TNNC1",
             "TNNI3", "TPM1", "TRIM63", "TNNT2", "ACTN2", "PLN", "FLNC", "PRKAG2", "CACNA1C", "LAMP2", "GLA"]


def sample_names(table):
    """Sample names are taken from the <sample>_gt columns of the extracted table."""
    header = pd.read_csv(table, sep="\t", nrows=0).columns
    return [c[:-3] for c in header if c.endswith("_gt")]


def labels(samples):
    """Study aliases A, B, C ... in column order."""
    return dict(zip(samples, string.ascii_uppercase))


def genotype(series):
    return series.str.replace("|", "/", regex=False).map(GT_CODE)


def allele_depths(series):
    ad = series.str.split(",", expand=True).iloc[:, :2].apply(pd.to_numeric, errors="coerce")
    total = ad[0] + ad[1]
    return ad[1], total


def load_sheet(workbook, name):
    """Rows of one worksheet as dictionaries (read-only, cached values)."""
    from openpyxl import load_workbook
    book = load_workbook(workbook, read_only=True, data_only=True)
    rows = iter(book[name].values)
    header = next(rows)
    data = [dict(zip(header, row)) for row in rows if any(v is not None for v in row)]
    book.close()
    return data
