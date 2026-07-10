#!/usr/bin/env python3

"""
mergeAdjascentAlns.py

Merge adjacent genome-alignment tiles in SAM/BAM/CRAM by inserting
query gaps as CIGAR insertions (I), using hard-clipping (H) as absolute
query coordinates.

This is intended for "genome-as-read" alignments where:
  - each alignment is a tile of a chromosome
  - leading H encodes query coordinate
  - alignments should be merged when adjacent on the reference

------------------------------------------------------------
Example

Input:

chr1_alt  0  Chr1  1000  60  100H10M5I5M20D50H
chr1_alt  0  Chr1  1035  60  135H20M3I10M40H

Step-by-step:

Alignment 1:
  query start = 100
  query span  = 10M + 5I + 5M = 20
  query end   = 120

Alignment 2:
  query start = 135

Insertion:
  gap = 135 - 120 = 15 bp

Merged CIGAR:
  100H10M5I5M20D15I20M3I10M40H

------------------------------------------------------------
Usage

  mergeAdjascentAlns.py [options] <input> <output>

Special values:
  - input  = "-" → read from stdin
  - output = "-" → write to stdout

------------------------------------------------------------
Examples

  mergeAdjascentAlns.py input.bam output.bam

  mergeAdjascentAlns.py input.bam - --output-fmt ubam \
      | samtools sort -o sorted.bam

  samtools view -h input.bam \
      | mergeAdjascentAlns.py - - \
      | samtools view -b > merged.bam

CRAM:

  mergeAdjascentAlns.py input.cram output.cram --reference ref.fa

------------------------------------------------------------
"""

import pysam
import re
import sys
from optparse import OptionParser

# ---------- CIGAR helpers ----------

cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')

def parse_cigar(cigar):
    return [(int(n), op) for n, op in cigar_re.findall(cigar)]

def query_span(cigar):
    return sum(n for n, op in parse_cigar(cigar)
               if op in ('M', 'I', '=', 'X'))

def leading_H(cigar):
    m = re.match(r'^(\d+)H', cigar)
    return int(m.group(1)) if m else 0

def strip_leading_H(cigar):
    return re.sub(r'^\d+H', '', cigar, count=1)

def strip_trailing_H(cigar):
    return re.sub(r'\d+H$', '', cigar, count=1)

# ---------- coordinate logic ----------

def query_coords(rec):
    cig = rec.cigarstring
    q_start = leading_H(cig)
    q_end = q_start + query_span(cig)
    return q_start, q_end

# ---------- merge condition ----------

def can_merge(a, b):
    if a.query_name != b.query_name:
        return False
    if a.reference_id != b.reference_id:
        return False
    if a.is_reverse != b.is_reverse:
        return False
    if a.reference_end != b.reference_start:
        return False

    a_qs, a_qe = query_coords(a)
    b_qs, _ = query_coords(b)

    if b_qs < a_qe:
        return False

    return True

# ---------- merging ----------

def merge_pair(a, b):
    a_qs, a_qe = query_coords(a)
    b_qs, _ = query_coords(b)

    gap = b_qs - a_qe
    if gap < 0:
        return None

    new = pysam.AlignedSegment(a.header)

    new.query_name = a.query_name
    new.reference_id = a.reference_id
    new.reference_start = a.reference_start
    new.flag = a.flag
    new.mapping_quality = min(a.mapping_quality, b.mapping_quality)

    new.cigarstring = (
        strip_trailing_H(a.cigarstring) +
        f"{gap}I" +
        strip_leading_H(b.cigarstring)
    )

    # optional sequence preservation
    new.query_sequence = a.query_sequence
    new.query_qualities = a.query_qualities

    return new

# ---------- PG header ----------

def add_pg_header(in_header, cmdline):
    header = in_header.to_dict()
    pg_lines = header.get("PG", [])

    base_id = "mergeAdjascentAlns"
    existing_ids = {pg["ID"] for pg in pg_lines} if pg_lines else set()

    pg_id = base_id
    i = 1
    while pg_id in existing_ids:
        pg_id = f"{base_id}.{i}"
        i += 1

    new_pg = {
        "ID": pg_id,
        "PN": "mergeAdjascentAlns",
        "VN": "1.0",
        "CL": cmdline
    }

    if pg_lines:
        new_pg["PP"] = pg_lines[-1]["ID"]

    header.setdefault("PG", []).append(new_pg)

    return pysam.AlignmentHeader.from_dict(header)

# ---------- stream merge ----------

def merge_stream(infile, outfile, debug=False):
    prev = None

    for rec in infile:
        if prev is None:
            prev = rec
            continue

        if can_merge(prev, rec):
            merged = merge_pair(prev, rec)
            if merged is not None:
                if debug:
                    a_qs, a_qe = query_coords(prev)
                    b_qs, _ = query_coords(rec)
                    sys.stderr.write(
                        f"MERGED {rec.query_name}: gap={b_qs - a_qe}\n")
                prev = merged
                continue

        if debug:
            sys.stderr.write(
                f"SKIP {rec.query_name}: not mergeable\n")

        outfile.write(prev)
        prev = rec

    if prev is not None:
        outfile.write(prev)

# ---------- I/O ----------

def open_alignment(path, mode, reference=None, template=None):
    if path == "-":
        handle = sys.stdin if "r" in mode else sys.stdout.buffer
        return pysam.AlignmentFile(
            handle,
            mode,
            reference_filename=reference,
            template=template
        )
    return pysam.AlignmentFile(
        path,
        mode,
        reference_filename=reference,
        template=template
    )

# ---------- main ----------

def main():
    parser = OptionParser(
        usage="%prog [options] <input> <output>",
        description="Merge adjacent alignment blocks using query-space gaps."
    )

    parser.add_option(
        "--reference", dest="reference",
        help="Reference FASTA (required for CRAM)"
    )

    parser.add_option(
        "--output-fmt", dest="outfmt", default="bam",
        help="Output format: sam | bam | ubam | cram (default: bam)"
    )

    parser.add_option(
        "--debug", action="store_true", dest="debug", default=False
    )

    options, args = parser.parse_args()

    if len(args) != 2:
        parser.error("Please provide input and output files")

    in_path, out_path = args

    # input mode
    in_mode = "r" if (in_path == "-" or in_path.endswith(".sam")) else "rb"

    # output mode
    fmt = options.outfmt.lower()
    if fmt == "sam":
        out_mode = "w"
    elif fmt == "bam":
        out_mode = "wb"
    elif fmt == "ubam":
        out_mode = "wbu"
    elif fmt == "cram":
        out_mode = "wc"
    else:
        parser.error("Unknown output format")

    infile = open_alignment(in_path, in_mode, options.reference)

    # add PG header
    cmdline = " ".join(sys.argv)
    new_header = add_pg_header(infile.header, cmdline)

    outfile = pysam.AlignmentFile(
        sys.stdout.buffer if out_path == "-" else out_path,
        out_mode,
        header=new_header,
        reference_filename=options.reference
    )

    merge_stream(infile, outfile, debug=options.debug)

    infile.close()
    outfile.close()


if __name__ == "__main__":
    main()
