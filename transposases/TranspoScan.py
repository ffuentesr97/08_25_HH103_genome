import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# Parameters
KEYWORDS = [
    "transposase",
    "insertion sequence",
    "mobile element",
    "IS family",
    "transposable element",
    "tnpA"
]

FLANK_SIZE = 150
MIN_IR_LEN = 12
MAX_IR_LEN = 35
IR_MAX_MISMATCH_RATIO = 0.15
TSD_MIN_LEN = 2
TSD_MAX_LEN = 10
TSD_MAX_MISMATCHES = 1
WINDOW_SIZE = 10000  # 10 kb

is_family_regex = re.compile(r'\bIS\d+\b', re.IGNORECASE)

def contains_keyword(text):
    text_lower = text.lower()
    if any(k.lower() in text_lower for k in KEYWORDS):
        return True
    if is_family_regex.search(text):
        return True
    return False

def count_mismatches(s1, s2):
    return sum(1 for a, b in zip(s1, s2) if a != b)

def detect_inverted_repeats(seq):
    seq = str(seq)
    for length in range(MAX_IR_LEN, MIN_IR_LEN - 1, -1):
        left = seq[:length]
        right = seq[-length:]
        right_rc = str(Seq(right).reverse_complement())
        mismatches = count_mismatches(left, right_rc)
        allowed_mismatches = max(1, int(length * IR_MAX_MISMATCH_RATIO))
        if mismatches <= allowed_mismatches:
            return f"Yes ({length} bp, {mismatches} mismatches)"
    return "No"

def detect_tsd(seq):
    seq = str(seq)
    for length in range(TSD_MAX_LEN, TSD_MIN_LEN - 1, -1):
        left = seq[:length]
        right = seq[-length:]
        mismatches = count_mismatches(left, right)
        if mismatches <= TSD_MAX_MISMATCHES:
            return f"Yes ({length} bp, {mismatches} mismatches)"
    return "Unclear"

def extract_transposases(gff_file, fasta_file,
                         out_fasta="transposases.fna",
                         out_tsv="transposase_info.tsv",
                         out_density="transposase_density_10kb.tsv",
                         out_summary="transposase_summary_per_replicon.tsv"):

    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    transposase_records = []
    info_table = []
    coordinates_by_replicon = defaultdict(list)

    with open(gff_file) as gff:
        for i, line in enumerate(gff):
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue

            seq_id, source, feature_type, start, end, score, strand, phase, attributes = parts
            if not contains_keyword(attributes):
                continue

            start, end = int(start), int(end)
            start_flank = max(0, start - FLANK_SIZE)
            end_flank = end + FLANK_SIZE

            sequence = sequences.get(seq_id)
            if not sequence:
                continue

            region = sequence.seq[start_flank:end_flank]
            if strand == "-":
                region = region.reverse_complement()

            record = SeqRecord(region, id=f"transposase_{i+1}", description=f"{seq_id}:{start}-{end}({strand})")
            transposase_records.append(record)

            ir_result = detect_inverted_repeats(region)
            tsd_result = detect_tsd(region)
            length_bp = abs(end - start) + 1
            annotation = attributes.split(";")[0]

            info_table.append({
                "ID": record.id,
                "Name": annotation,
                "Start": start,
                "End": end,
                "Strand": strand,
                "Length (bp)": length_bp,
                "Inverted Repeats": ir_result,
                "TSDs": tsd_result,
                "Annotation": annotation
            })

            coordinates_by_replicon[seq_id].append((start, end))

    # Write transposase sequences to FASTA
    SeqIO.write(transposase_records, out_fasta, "fasta")

    # Write detailed TSV table
    with open(out_tsv, "w") as out:
        headers = ["ID", "Name", "Start", "End", "Strand", "Length (bp)", "Inverted Repeats", "TSDs", "Annotation"]
        out.write("\t".join(headers) + "\n")
        for row in info_table:
            out.write("\t".join(str(row[col]) for col in headers) + "\n")

    # Compute transposase density per 10 kb window
    with open(out_density, "w") as out:
        out.write("Replicon\tStart\tEnd\tTransposases\n")
        for replicon, sequence in sequences.items():
            seq_length = len(sequence)
            hits = coordinates_by_replicon.get(replicon, [])
            for window_start in range(0, seq_length, WINDOW_SIZE):
                window_end = min(window_start + WINDOW_SIZE, seq_length)
                count = sum(1 for s, e in hits if not (e < window_start or s > window_end))
                out.write(f"{replicon}\t{window_start}\t{window_end}\t{count}\n")

    # Summary per replicon
    with open(out_summary, "w") as summary:
        summary.write("Replicon\tTotal_Transposases\tReplicon_Length_bp\tProportion_per_10kb\n")
        for replicon, sequence in sequences.items():
            replicon_length = len(sequence)
            total_hits = len(coordinates_by_replicon.get(replicon, []))
            proportion = round((total_hits / replicon_length) * 10000, 3) if replicon_length > 0 else 0
            summary.write(f"{replicon}\t{total_hits}\t{replicon_length}\t{proportion}\n")

    # Final output messages
    print(f"Transposases detected: {len(transposase_records)}")
    print(f"FASTA output: {out_fasta}")
    print(f"Structural details: {out_tsv}")
    print(f"Density per 10 kb: {out_density}")
    print(f"Summary per replicon: {out_summary}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python detect_transposases.py input.gff genome.fna")
        sys.exit(1)
    gff_input, fasta_input = sys.argv[1], sys.argv[2]
    extract_transposases(gff_input, fasta_input)
