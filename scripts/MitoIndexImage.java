import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

/** Build a GATK BWA image from verified existing indexes, without reindexing FASTA. */
public final class MitoIndexImage {
    public static void main(String[] args) throws Exception {
        if (args.length != 3) throw new IllegalArgumentException("indexPrefix output.img reference.fai");
        if (Files.exists(Paths.get(args[1]))) throw new IllegalArgumentException("Output already exists");
        BwaMemIndex.createIndexImageFromIndexFiles(args[0], args[1]);
        List<String> expected = new ArrayList<>();
        for (String line : Files.readAllLines(Paths.get(args[2]))) expected.add(line.split("\t")[0]);
        try (BwaMemIndex image = new BwaMemIndex(args[1])) {
            if (!image.getReferenceContigNames().equals(expected)) {
                throw new IllegalStateException("Image/reference sequence dictionary mismatch");
            }
            System.out.println("{\"status\":\"PASS\",\"contig_names_match\":true,\"contigs\":"
                    + expected.size() + ",\"native_bwa_version\":\"" + BwaMemIndex.getBWAVersion() + "\"}");
        }
    }
}
