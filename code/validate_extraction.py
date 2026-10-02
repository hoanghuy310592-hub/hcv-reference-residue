#!/usr/bin/env python3
"""Independent re-extraction of the reference residues, as a check on the
primary BioEdit analysis described in Supplementary Methods S3.

Method: translate each genome, anchor the polyprotein frame on the NS5B GDD
catalytic motif (NS5B residue 317), align each polyprotein to the subtype 1a
reference with a BLOSUM62 global alignment, and read the fourteen positions
through that alignment. No multiple-sequence alignment is used, so the two
analyses share no step beyond the input sequences.

Usage:  python code/validate_extraction.py
"""
import csv, json, re, sys
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices

FASTA   = "data/hcv_93_subtype_references.fas"
MATRIX  = "data/fig2_full_residue_matrix.csv"
OUT     = "data/validation_independent_extraction.csv"

AA = "ACDEFGHIKLMNPQRSTVWYX"
POS_NS5A = [24, 28, 30, 31, 32, 38, 58, 62, 92, 93]
POS_NS5B = [159, 282, 316, 321]
NS5A_START, NS5B_START = 1973, 2421        # H77 polyprotein numbering
GDD_POLY = NS5B_START + 316                # GDD at NS5B residue 317
REF_ACC  = "M62321.1"


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip().upper()
    return seqs


CORE = {"GGGQIVGG": 26, "GGVYLLPR": 33, "GVYLLPR": 34, "QIVGG": 29}


def _trim(prot, start, minlen):
    if start < 0:
        return None
    tail = prot[start:]
    stop = tail.find("*")
    if stop != -1:
        tail = tail[:stop]
    return tail if len(tail) >= minlen else None


def polyprotein(nt):
    """Anchor the polyprotein frame on the NS5B GDD motif, falling back to a
    conserved core motif for the few genomes whose GDD search is ambiguous."""
    best = None
    for frame in range(3):
        sub = nt[frame:len(nt) - ((len(nt) - frame) % 3)]
        prot = "".join(c if c in AA else "X" for c in str(Seq(sub).translate()))
        for m in re.finditer("GDD", prot):
            cand = _trim(prot, m.start() - (GDD_POLY - 1), 2990)
            if cand and (best is None or len(cand) > len(best)):
                best = cand
    if best is not None:
        return best
    for frame in range(3):
        sub = nt[frame:len(nt) - ((len(nt) - frame) % 3)]
        prot = "".join(c if c in AA else "X" for c in str(Seq(sub).translate()))
        for mot, off in CORE.items():
            i = prot.find(mot)
            if i < 0:
                continue
            cand = _trim(prot, i - (off - 1), 2900)
            if cand and (best is None or len(cand) > len(best)):
                best = cand
    return best


def main():
    seqs = read_fasta(FASTA)
    ref = polyprotein(seqs[REF_ACC])
    if ref is None:
        sys.exit("could not anchor the reference polyprotein")

    aligner = Align.PairwiseAligner(
        mode="global", open_gap_score=-11, extend_gap_score=-1,
        substitution_matrix=substitution_matrices.load("BLOSUM62"))

    def mapping(a_seq, b_seq):
        aln = aligner.align(a_seq, b_seq)[0]
        out, i, j = {}, 0, 0
        for x, y in zip(aln[0], aln[1]):
            if x != "-" and y != "-":
                out[i] = y
            if x != "-":
                i += 1
            if y != "-":
                j += 1
        return out

    cols = [f"NS5A-{p}" for p in POS_NS5A] + [f"NS5B-{p}" for p in POS_NS5B]
    rows, failed = [], []
    for acc, nt in seqs.items():
        prot = polyprotein(nt)
        if prot is None:
            failed.append(acc)
            continue
        m = mapping(ref, prot)
        rec = {"accession": acc, "polyprotein_length": len(prot),
               "gdd_at_ns5b_residue": 317}
        for p in POS_NS5A:
            rec[f"NS5A-{p}"] = m.get(NS5A_START - 2 + p, "?")
        for p in POS_NS5B:
            rec[f"NS5B-{p}"] = m.get(NS5B_START - 2 + p, "?")
        rows.append(rec)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    primary = {r["accession"]: r for r in csv.DictReader(open(MATRIX))}
    by_acc = {r["accession"]: r for r in rows}
    total = agree = 0
    diffs = []
    for acc, p in primary.items():
        q = by_acc.get(acc)
        if q is None:
            continue
        for c in cols:
            total += 1
            if p[c] == q[c]:
                agree += 1
            else:
                diffs.append((p["subtype"], acc, c, p[c], q[c]))

    print(f"genomes anchored        {len(rows)} / {len(seqs)}")
    if failed:
        print(f"not anchored            {', '.join(failed)}")
    print(f"cells compared          {total}")
    print(f"agreement               {agree} ({100 * agree / total:.2f}%)")
    print(f"disagreement            {len(diffs)}")
    for d in diffs:
        print("   subtype %-5s %-12s %-10s primary=%s  independent=%s" % d)


if __name__ == "__main__":
    main()
