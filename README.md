Multi-CRAM utils
================

This contains various scripts written by Copilot upon my request, to
post-process the CRAM files produced by the
[`nf-core/pairgenomealign`](https://nf-co.re/pairgenomealign/) pipeline.

The CRAM files represent pairwise whole-genome alignments of mutliple _query_
genomes to a single _target_ genome.  The _target_ genome takes the role of the
reference sequence in the CRAM files.  The _query_ genomes are represented as
read groups.  Hard-clipped regions can be longer than the 268,435,455 length
limit (28 bits) in the BAM and CRAM files; in CIGAR lines they are represented
as adjascent clipping operations.

The `mergeAdjascentAlns.sam` is toy test data.
