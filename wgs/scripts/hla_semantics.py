"""Shared HLA allele-name and KIR-ligand semantics for steps 24 and 25. Standard library only."""
import re

# Bw4/Bw6 follow the KIR-ligand families of Abraham et al. (PLoS Genet 2008): the B*13/27/37/38/
# 44/47/49/51/52/53/57/58/59/63/77 groups are Bw4 throughout. B*15 is split (15:10 and a few
# others are Bw4, most are Bw6), so an unmatched B*15 allele is reported unknown, not guessed.
BW4_FAMILIES = {"13", "27", "37", "38", "44", "47", "49", "51", "52", "53", "57", "58", "59", "63", "77"}
BW4_B15 = {"15:10"}


def two_fields(allele):
    """Keep the first two colon fields of an HLA allele name (02:07:01 -> 02:07).

    rsplit(':', 1)[0] truncated a two-field allele to one field (44:02 -> 44), after which no
    two-field whitelist entry could ever match. Returns the input unchanged when it has no
    digit:digit prefix.
    """
    m = re.match(r"(\d+:\d+)", str(allele))
    return m.group(1) if m else str(allele)


def bw_epitope(allele_two_field):
    """Bw4 / Bw6 / unknown for a two-field B allele, by KIR-ligand family."""
    fam = allele_two_field.split(":")[0]
    if fam == "15":
        return "Bw4" if allele_two_field in BW4_B15 else "unknown (B*15 split family)"
    return "Bw4" if fam in BW4_FAMILIES else "Bw6"


def c_ligand(allele_two_field):
    """C1 / C2 / unknown for a two-field C allele (Lys80 = C1, Asn80 = C2 groupings)."""
    fam = allele_two_field.split(":")[0]
    c1 = {"01", "03", "07", "08", "12", "14", "16"}
    c2 = {"02", "04", "05", "06", "15", "17", "18"}
    return "C1" if fam in c1 else "C2" if fam in c2 else "?"


def a311(alleles_two_field):
    """The HLA-A alleles among KIR3DL2's A3/A11 ligand groups."""
    return [a for a in alleles_two_field if a.split(":")[0] in ("03", "11")]
