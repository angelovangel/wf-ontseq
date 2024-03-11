#!/usr/bin/env python3

# view genbank record with coverage

import sys
import os
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from Bio import SeqIO
#import numpy as np
import pandas as pd

# custom coloring 
class MyTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            return  "#A9CCE3"
        elif feature.type == "promoter":
            return "#1E8449"
        elif feature.type == "terminator":
            return "#F0B27A"
        elif feature.type == "misc_feature":
            return "#B2BABB"
        else:
            return "#F1948A"


# argv[1] is gbk, argv[2] is coverage from perbase
gbk = sys.argv[1]
samplename = os.path.basename(gbk).split('.', 1)[0]
covfile = sys.argv[2]
outfile = os.path.dirname(covfile) + '/' + samplename + '.coverage.svg'


fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(12, 6), sharex=True, gridspec_kw={"height_ratios": [5, 1]}
)

# read gbk and plot
record = SeqIO.read(gbk, format="genbank")
graphic_record = MyTranslator().translate_record(record)
graphic_record.plot(ax=ax1, with_ruler=True, strand_in_label_threshold=10)

# read cov and plot
tsv = pd.read_csv(covfile, sep="\t", index_col=0)
cov = tsv["DEPTH"]
#print(cov.values)

ax2.fill_between(x = tsv["POS"], y1 = 0, y2 = cov.values, alpha=0.3)
#ax2.plot(tsv["POS"], cov.values, linewidth = 1)
ax2.set_ylim(bottom=0)
ax2.set_ylabel("Coverage")
ax2.spines['top'].set_visible(False)
fig.suptitle("Features and read coverage for " + str(samplename))
fig.tight_layout()
fig.savefig(outfile)
