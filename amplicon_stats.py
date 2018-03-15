#!/usr/bin/env python

import argparse
from collection import defaultdict
import matplotlib
import re


def plot_depth(bams, locus):
    """Plots depth across given loci."""
    plot_map = defaultdict(list)
    for bam in bams:
        samplename = bam[re.search(".bam", bam).start() : ]
        chrom, begin, end = locus
        


def main():
    """"""
    # arg parsing
    desc = ""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--loci", required=True,
                        nargs=3, metavar=("chr", "begin", "end"),
                        action="append", help="comma-delimited loci [required]")
    parser.add_argument("-b", "--bam", metavar="BAM", nargs='+',
                        required=True, help="BAM file(s) [required]")
    args = parser.parse_args()
    # central dispatch
    #plot_depth()


if __name__ == "__main__":
    main()
