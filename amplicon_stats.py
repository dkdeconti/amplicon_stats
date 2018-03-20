#!/usr/bin/env python

import argparse
from collection import defaultdict
import pysam
import matplotlib
import re


def plot_depth(pileups, locus):
    """Plots depth across given loci."""
    plot_map = defaultdict(list)
    #
        
def parse_bam_for_pileup(bams, loci):
    """Parses BAM for pileup data into map by locus key."""
    pileups = {}
    for locus in loci:
        chrom = locus[0]
        begin = int(locus[1])
        end = int(locus[2])
        pileups[loci] = [get_pileup_vector(bam, chrom, begin, end)
                         for bam in bams]
    return pileups


def get_pileup_vector(bam, chrom, begin, end):
    """Get depth vector from bam given pos."""
    return [p.nsegments
            for p in pysam.AlignmentFile(bam, 'rb').pileup(chrom, begin, end)]


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
