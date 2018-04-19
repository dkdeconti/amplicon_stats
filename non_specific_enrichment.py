#! /usr/bin/env python

import argparse
import pysam
import numpy
import random
from scipy.stats.mstats import gmean


def get_mean_depth(bam, ref, start, end):
    """
    """
    depths = [0]*(end - start)
    for pileupcolumn in bam.pileup(ref, start, end + 1):
        if pileupcolumn.reference_pos >= end:
            continue
        depths[pileupcolumn.reference_pos - start] = pileupcolumn.nsegments
    #return gmean(depths)
    return numpy.mean(depths)
    #return numpy.median(depths)


def get_stats(bam, bin_size):
    """
    Calculates geometric mean of depth over tiling windows across genome to
    return mean and standard deviation.

    params
    bamfilename: filename (and path) for BAM file (string)
    """
    mean_depths = [] # init the first
    refs = zip(bam.references, bam.lengths)
    for i in range(100):
        #if i % 1000 == 0:
        #    print i
        ref_len = 0
        while ref_len < bin_size:
            ref, ref_len = random.choice(refs)
        start = random.randint(0, ref_len - bin_size)
        mean_depths.append(get_mean_depth(bam, ref, start, start + bin_size))
    #print mean_depths
    print numpy.nanmean(mean_depths)


def main():
    """
    """
    # arg parsing
    desc = ""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-l", "--loci", required=True,
                        nargs=3, metavar=("chr", "begin", "end"),
                        help="comma-delimited loci [required]")
    parser.add_argument("-b", "--bam", metavar="BAM",
                        required=True, help="BAM file [required]")
    args = parser.parse_args()
    # central dispatch
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        bin_size = int(args.loci[2]) - int(args.loci[1])
        ref, start, end = args.loci
        start, end = int(start), int(end)
        print get_mean_depth(bam, ref, start, end)
        get_stats(bam, bin_size)


if __name__ == "__main__":
    main()
