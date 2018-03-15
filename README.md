# amplicon_stats
From "cut" DNA, identify specificity of enrichment of "cut" amplicon over background, partuclarly around cut loci.

## Summary

Creates an HTML report of enrichment of amplicon reads by a [redacted]-based
cutting method for genomic DNA. Loci of cutting are specific to the base
pair level. By providing the loci of the cut sites and BAM files, this
creates an HTML report with:
* Depth plot across the loci +/- a given BP range
* Depth plot at the locus +/- a given BP range
* Identified potential alternative sites of (non-specific) enriched regions

## Usage
