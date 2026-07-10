#!/usr/bin/env python3
import pysam
import re
import sys
import signal

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

cigar_re = re.compile(r"(\d+)([MIDNSHP=X])")
md_re = re.compile(r"(\d+|\^[A-Z]+|[A-Z])")


def parse_cigar(cigar):
    return [(int(n), op) for n, op in cigar_re.findall(cigar)]


def iterate_md(md):
    return md_re.findall(md)


def extract_records(rec, fasta):
    ref_seq = fasta.fetch(rec.reference_name)
    cigar = parse_cigar(rec.cigarstring)
    md_tokens = iterate_md(rec.get_tag("MD")) if rec.has_tag("MD") else None

    ref_pos = rec.reference_start
    read_pos = 0
    seq = rec.query_sequence

    md_i = 0
    md_remaining = 0

    def next_md():
        nonlocal md_i, md_remaining
        while md_i < len(md_tokens):
            tok = md_tokens[md_i]
            md_i += 1
            if tok.isdigit():
                md_remaining = int(tok)
                return "match"
            elif tok.startswith("^"):
                return ("deletion", tok[1:])
            else:
                return ("mismatch", tok)
        return None

    md_state = next_md() if md_tokens else None

    # per-position aggregation
    pos_buffer = {}

    def record(pos, ref, alt, is_variant):
        if pos not in pos_buffer:
            pos_buffer[pos] = {"ref": ref, "variant": None}
        else:
            # ✅ Do NOT overwrite an existing variant
            if pos_buffer[pos]["variant"] is not None:
                return

            # ✅ Keep consistent REF (avoid corruption)
            if len(ref) == 1 and pos_buffer[pos]["ref"] != ref:
                return

        if is_variant:
            if alt != ref:
                pos_buffer[pos]["variant"] = alt

    for length, op in cigar:
        if op in ("M", "=", "X"):
            for _ in range(length):
                if md_state == "match":
                    if md_remaining > 0:
                        pos = ref_pos + 1
                        # ✅ only emit match if no variant already recorded
                        if pos not in pos_buffer:
                            ref_base = ref_seq[ref_pos]
                            record(pos, ref_base, "<*>", False)
                        md_remaining -= 1
                        ref_pos += 1
                        read_pos += 1
                        continue
                    else:
                        md_state = next_md()

                if isinstance(md_state, tuple) and md_state[0] == "mismatch":
                    ref_base = ref_seq[ref_pos]
                    alt_base = seq[read_pos]

                    if alt_base == ref_base:
                        ref_pos += 1
                        read_pos += 1
                        md_state = next_md()
                        continue

                    record(ref_pos + 1, ref_base, alt_base, True)

                    ref_pos += 1
                    read_pos += 1
                    md_state = next_md()

                elif isinstance(md_state, tuple) and md_state[0] == "deletion":
                    # deletion handled via CIGAR
                    md_state = next_md()

                else:
                    pos = ref_pos + 1
                    if pos not in pos_buffer:
                        ref_base = ref_seq[ref_pos]
                        record(pos, ref_base, "<*>", False)

                    ref_pos += 1
                    read_pos += 1

        elif op == "I":
            alt_seq = seq[read_pos:read_pos + length]
            anchor = ref_seq[ref_pos - 1]

            record(ref_pos, anchor, anchor + alt_seq, True)

            read_pos += length

        elif op == "D":
            deleted = ref_seq[ref_pos:ref_pos + length]
            anchor = ref_seq[ref_pos - 1]

            # ✅ correct VCF deletion
            record(ref_pos, anchor + deleted, anchor, True)

            ref_pos += length

        elif op in ("S", "H"):
            continue

        elif op in ("N", "P"):
            ref_pos += length

    # emit gVCF-style
    for pos in sorted(pos_buffer):
        ref = pos_buffer[pos]["ref"]
        var = pos_buffer[pos]["variant"]

        if var is not None:
            yield rec.reference_name, pos, ref, var
        else:
            yield rec.reference_name, pos, ref, "<*>"


def write_header(out, fasta):
    out.write("##fileformat=VCFv4.2\n")
    for c in fasta.references:
        out.write(f"##contig=<ID={c},length={fasta.get_reference_length(c)}>\n")
    out.write("##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n")
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def main():
    bam = sys.argv[1]
    ref = sys.argv[2]

    if bam == "-":
        aln = pysam.AlignmentFile(sys.stdin)
    else:
        aln = pysam.AlignmentFile(bam)

    fasta = pysam.FastaFile(ref)

    write_header(sys.stdout, fasta)

    for rec in aln:
        if rec.is_unmapped:
            continue

        try:
            for chrom, pos, ref_base, alt_base in extract_records(rec, fasta):
                sys.stdout.write(
                    f"{chrom}\t{pos}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n"
                )
        except Exception:
            continue


if __name__ == "__main__":
    main()
