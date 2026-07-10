#!/usr/bin/env python3

import argparse
import re
import sys

import pysam


def parse_region(region: str):
    """
    Parse region strings like:
      scaffold_1:10000000-10000100
      scaffold_1:10,000,000-10,000,100

    Returns:
      chrom, start_1based, end_1based
    """
    m = re.fullmatch(r"([^:]+):([\d,]+)-([\d,]+)", region)
    if not m:
        raise ValueError(
            f"Invalid region format: {region!r}. "
            "Expected something like scaffold_1:10000000-10000100"
        )

    chrom = m.group(1)
    start = int(m.group(2).replace(",", ""))
    end = int(m.group(3).replace(",", ""))

    if start < 1:
        raise ValueError("Region start must be >= 1")
    if end < start:
        raise ValueError("Region end must be >= start")

    return chrom, start, end


def wrapped(seq: str, width: int = 80):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract only the query bases aligned inside a genomic region "
            "from a CRAM file, and print them as FASTA."
        )
    )
    parser.add_argument("--cram", required=True, help="Input CRAM file")
    parser.add_argument("--ref", required=True, help="Reference FASTA used by the CRAM")
    parser.add_argument(
        "--region",
        required=True,
        help="Region in the form scaffold_1:10000000-10000100",
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output FASTA file (default: stdout)",
    )
    parser.add_argument(
        "--include-deletions",
        action="store_true",
        help=(
            "Include deletions relative to the reference as '-' in the output. "
            "Default: deletions are skipped."
        ),
    )
    args = parser.parse_args()

    try:
        chrom, start_1based, end_1based = parse_region(args.region)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")

    # pysam fetch uses 0-based, half-open coordinates
    fetch_start = start_1based - 1
    fetch_end = end_1based

    try:
        aln_file = pysam.AlignmentFile(
            args.cram,
            "rc",
            reference_filename=args.ref,
        )
    except Exception as e:
        sys.exit(f"ERROR: could not open CRAM/reference: {e}")

    out = sys.stdout if args.output == "-" else open(args.output, "w")

    n_written = 0

    try:
        for aln in aln_file.fetch(chrom, fetch_start, fetch_end):
            if aln.is_unmapped:
                continue
            if aln.query_sequence is None:
                continue

            seq = []
            first_ref = None
            last_ref = None

            # get_aligned_pairs gives 0-based reference positions
            for qpos, rpos in aln.get_aligned_pairs(matches_only=False):
                if rpos is None:
                    continue

                if fetch_start <= rpos < fetch_end:
                    if first_ref is None:
                        first_ref = rpos + 1
                    last_ref = rpos + 1

                    if qpos is None:
                        if args.include_deletions:
                            seq.append("-")
                    else:
                        seq.append(aln.query_sequence[qpos])

            if seq:
                try:
                    rg = aln.get_tag("RG")
                except KeyError:
                    rg = "RG:NA"

                description = (
                    f"query_name={aln.query_name}"
                    f"|region={chrom}:{first_ref}-{last_ref}"
                    f"|flag={aln.flag}"
                    f"|mapq={aln.mapping_quality}"
                )

                print(f">{rg} {description}", file=out)
                for line in wrapped("".join(seq)):
                    print(line, file=out)
                n_written += 1

    except ValueError as e:
        sys.exit(f"ERROR while fetching region {args.region}: {e}")
    finally:
        aln_file.close()
        if out is not sys.stdout:
            out.close()

    print(f"[INFO] Wrote {n_written} FASTA entries", file=sys.stderr)


if __name__ == "__main__":
    main()
