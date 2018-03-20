#!/usr/bin/env python

import argparse
from collection import defaultdict
import numpy
import matplotlib.pyplot as plt
import pysam
import matplotlib
import re


def plot_depth(pileups, locus, bams):
    """Plots depth across given loci."""
    chrom, begin, end = locus
    begin = int(begin)
    plot_map = defaultdict(list)
    filenames = [re.sub('\.bam$', '', bam) for bam in bams]
    plot_file = '-'.join(locus) + '.png'
    for i, pileup in enumerate(pileups):
        # TODO add a filepath
        fig, axarray = plt.subplots(nrows=len(pileups),
                                    ncols=1,
                                    figsize=(8, 2*len(pileups)),
                                    sharex=True)
        axarray[i].plot(range(begin, len(pileups) + begin),
                        pileup)
        axarray[i].set_title(filenames[i])
    fig.savefig(plot_file)
    fig.close()

        
def parse_bam_for_pileup(bams, loci, buf):
    """Parses BAM for pileup data into map by locus key."""
    pileups = {}
    for locus in loci:
        chrom = locus[0]
        begin = int(locus[1]) -1
        end = int(locus[2])
        pileups[loci] = [get_pileup_vector(bam, chrom, begin, end, buf)
                         for bam in bams]
    return pileups


def get_pileup_vector(bam, chrom, begin, end, buf):
    """Get depth vector from bam given pos."""
    return [p.nsegments
            for p in pysam.AlignmentFile(bam, 'rb').pileup(chrom,
                                                           begin - buf,
                                                           end + buf)]


def main():
    """"""
    # arg parsing
    desc = ""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--loci", required=True,
                        nargs=3, metavar=("chr", "begin", "end"),
                        action="append", help="comma-delimited loci [required]")
    parser.add_argument("--buf", metavar="BUFFER", type=int, default=1000,
                        help="buffer about the loci to include in plots")
    parser.add_argument("-b", "--bam", metavar="BAM", nargs='+',
                        required=True, help="BAM file(s) [required]")
    args = parser.parse_args()
    # central dispatch
    parse_bam_for_pileup(args.bams, args.loci, args.buf)
    #plot_depth()


if __name__ == "__main__":
    main()