"""Shared VCF GT parser for the panel steps (11_pgx_extra, 27_candidate_genes). Standard library only."""


def gt_alleles(gt, ref, alt_field):
    """Split a VCF GT into per-allele names against REF and the record's comma-separated ALT list.

    Returns a list with one entry per allele in the GT: the base/allele string, or None for a
    missing ('.') allele. Returns None for a full no-call ('./.', '.|.' or '.'). Allele indices
    beyond the record's ALT list render as '<altN>' instead of being silently mislabeled, and a
    token that is not a VCF allele index surfaces the raw GT verbatim. Both keep uninterpretable
    genotypes visible rather than letting them collapse into a reference homozygote, which is
    what a simple count of '1' characters did.
    """
    if gt in (".", "./.", ".|."):
        return None
    sep = "|" if "|" in gt else "/"
    alts = alt_field.split(",")
    out = []
    for tok in gt.split(sep):
        if tok == ".":
            out.append(None)
            continue
        try:
            i = int(tok)
        except ValueError:
            return ["GT:" + gt]
        if i == 0:
            out.append(ref)
        elif i <= len(alts):
            out.append(alts[i - 1])
        else:
            out.append(f"<alt{i}>")
    return None if all(a is None for a in out) else out
