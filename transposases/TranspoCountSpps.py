#!/usr/bin/env python3

import os
import re
import sys

# ------------------------------------------------------------------
# Keyword / pattern configuration for Prokka annotations
# ------------------------------------------------------------------
KEYWORDS = [
    r"transposase",
    r"insertion sequence",
    r"mobile[ _]?element",
    r"mobile[ _]?element[ _]?protein",
    r"IS\d+",
    r"tnpA",
]

PATTERN = re.compile(
    r"\b(" + r"|".join(KEYWORDS) + r")\b",
    re.IGNORECASE
)

# ------------------------------------------------------------------
# Count lines matching keywords in attribute field (column 9 of GFF)
# ------------------------------------------------------------------
def count_transposases_in_gff(gff_path):
    count = 0
    with open(gff_path) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attributes = fields[8]
            if PATTERN.search(attributes):
                count += 1
    return count

# ------------------------------------------------------------------
# Prefer GCF assemblies over GCA when both are present
# ------------------------------------------------------------------
def filter_prefer_gcf(results):
    assemblies = {}
    for name, total in results:
        match = re.search(r"(GC[AF]_\d{6,9})\.\d+", name)
        if not match:
            assemblies[name] = (name, total)
            continue
        base_accession = match.group(1)
        if base_accession.startswith("GCF_"):
            assemblies[base_accession] = (name, total)
        else:  # GCA
            gcf_version = base_accession.replace("GCA_", "GCF_")
            if gcf_version not in assemblies:
                assemblies[base_accession] = (name, total)
    return list(assemblies.values())

# ------------------------------------------------------------------
# Recursively find all GFF/GFF3 files in a directory
# ------------------------------------------------------------------
def find_gff_files_recursive(folder_path):
    gff_files = []
    for root, dirs, files in os.walk(folder_path):
        for fn in files:
            if fn.lower().endswith((".gff", ".gff3")):
                full_path = os.path.join(root, fn)
                gff_files.append(full_path)
    return sorted(gff_files)

# ------------------------------------------------------------------
# Analyze all GFFs in a directory, count hits, and write TSV summary
# ------------------------------------------------------------------
def analyze_directory(folder_path, output_file="transposase_counts.tsv"):
    gff_files = find_gff_files_recursive(folder_path)
    results = []
    for path in gff_files:
        relative_name = os.path.relpath(path, folder_path)
        total = count_transposases_in_gff(path)
        results.append((relative_name, total))
    filtered_results = filter_prefer_gcf(results)
    with open(output_file, "w") as out:
        out.write("Bacterium\tTransposase_or_MobileElement_Count\n")
        for name, total in filtered_results:
            out.write(f"{name}\t{total}\n")
    print(f"Done: {output_file}")

# ------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_transposases_gffs.py <gff_folder>")
        sys.exit(1)
    folder = sys.argv[1]
    analyze_directory(folder)

