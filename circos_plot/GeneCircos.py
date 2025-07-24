import sys
import os
import numpy as np
from Bio import SeqIO
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pycirclize import Circos
from pycirclize.parser import Genbank, Gff

# ----------------------------
# Validate command-line arguments
# ----------------------------
if len(sys.argv) != 4:
    print("Usage: python plot_mobile_elements.py <genbank_file.gbk> <genome_name> <annotations.gff>")
    sys.exit(1)

gbk_file = sys.argv[1]
genome_name = sys.argv[2]
gff_file = sys.argv[3]

# ----------------------------
# Load GenBank and GFF files
# ----------------------------
gbk = Genbank(gbk_file)
gff = Gff(gff_file)

# ----------------------------
# Keywords to identify mobile elements
# ----------------------------
keywords = [
    "transposase",
    "insertion sequence",
    "mobile element",
    "IS family",
    "transposable element",
    "tnpA"
]

def contains_mobile_keyword(qualifiers):
    """Returns True if any mobile-related keyword is found in feature annotations."""
    for key, values in qualifiers.items():
        if isinstance(values, str):
            values = [values]
        for value in values:
            if any(kw.lower() in value.lower() for kw in keywords):
                return True
    return False

# ----------------------------
# Create Circos plot and genome sector
# ----------------------------
genome_size_bp = sum(len(record) for record in gbk.records)
circos = Circos(sectors={gbk.name: gbk.range_size})
circos.text(f"$\\it{{Sinorhizobium\\ fredii}}$ HH103\n{genome_name}\n{genome_size_bp} bp", size=12, r=16)
sector = circos.get_sector(gbk.name)

# ----------------------------
# Outer tick track (scale)
# ----------------------------
outer_track = sector.add_track((98, 99))
outer_track.axis(fc="dimgrey")
outer_track.xticks_by_interval(50000, label_formatter=lambda v: f"{int(v / 1_000)} Kb")
outer_track.xticks_by_interval(25000, tick_length=1, show_label=False)

# ----------------------------
# CDS tracks (forward and reverse)
# ----------------------------
fwd_cds_track = sector.add_track((90, 97), r_pad_ratio=0.1)
fwd_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), fc="tomato")

rev_cds_track = sector.add_track((83, 90), r_pad_ratio=0.1)
rev_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), fc="deepskyblue")

# ----------------------------
# Highlight mobile genetic elements (CDS with keywords)
# ----------------------------
highlight_fwd_track = sector.add_track((74, 79), r_pad_ratio=0.1)
for feature in gff.extract_features("CDS", target_strand=1):
    if contains_mobile_keyword(feature.qualifiers):
        highlight_fwd_track.genomic_features(feature, fc="orange", ec="darkorange", lw=0.5)

highlight_rev_track = sector.add_track((74, 79), r_pad_ratio=0.1)
for feature in gff.extract_features("CDS", target_strand=-1):
    if contains_mobile_keyword(feature.qualifiers):
        highlight_rev_track.genomic_features(feature, fc="orange", ec="darkorange", lw=0.5)

# ----------------------------
# GC Content track
# ----------------------------
gc_content_track = sector.add_track((55, 70))
pos_list, gc_contents = gbk.calc_gc_content()
gc_contents = gc_contents - gbk.calc_genome_gc_content()
positive_gc = np.where(gc_contents > 0, gc_contents, 0)
negative_gc = np.where(gc_contents < 0, gc_contents, 0)
abs_gc = np.max(np.abs(gc_contents))
gc_content_track.fill_between(pos_list, positive_gc, 0, vmin=-abs_gc, vmax=abs_gc, color="black")
gc_content_track.fill_between(pos_list, negative_gc, 0, vmin=-abs_gc, vmax=abs_gc, color="grey")

# ----------------------------
# GC Skew track
# ----------------------------
gc_skew_track = sector.add_track((40, 55))
pos_list, gc_skews = gbk.calc_gc_skew()
positive_skew = np.where(gc_skews > 0, gc_skews, 0)
negative_skew = np.where(gc_skews < 0, gc_skews, 0)
abs_skew = np.max(np.abs(gc_skews))
gc_skew_track.fill_between(pos_list, positive_skew, 0, vmin=-abs_skew, vmax=abs_skew, color="olive")
gc_skew_track.fill_between(pos_list, negative_skew, 0, vmin=-abs_skew, vmax=abs_skew, color="purple")

# ----------------------------
# Plot and legend
# ----------------------------
fig = circos.plotfig()

legend_elements = [
    Patch(color="tomato", label="Forward CDS"),
    Patch(color="deepskyblue", label="Reverse CDS"),
    Patch(color="orange", label="Mobile elements"),
    Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
    Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None"),
    Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=6, ls="None"),
    Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None"),
]
_ = circos.ax.legend(handles=legend_elements, bbox_to_anchor=(0.52, 0.44), loc="center", fontsize=8)

# ----------------------------
# Save figure
# ----------------------------
output_img = f"result_plot_{genome_name}.png"
fig.savefig(output_img, dpi=1200)
print(f"Image saved: {output_img}")
