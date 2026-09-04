#!/usr/bin/env python3
"""
design_3primer_genotyping.py

Version: v1
Author: Sim Sakong, Hansen Lab, MIT
Date: 2026-09-04

Automated 3-primer design for CRISPR genotyping by PCR product size — covers
both a deletion/knockout (KO) and an HDR knock-in (KI) edit, since the two
turn out to be the exact mirror image of each other (see the "Knock-in"
section further down).

*** This script does NOT implement its own Tm/GC/complementarity calculations. ***
It calls the installed primer3-py package (Python bindings to the real Primer3 C
engine) directly to compute Tm, GC%, self-complementarity, hairpin, and primer-primer
dimer scores. This is the same engine (SantaLucia nearest-neighbor thermodynamics +
the thal algorithm) used internally by Primer3Plus / Primer-BLAST.

Setup (one time only, inside your venv):
    pip install primer3-py

primer3-py bundles its own thermodynamic parameters, so there's no need to compile
primer3_core or ntthal binaries separately.

Core idea (3-primer genotyping), in generic terms:
    Every edit (deletion or insertion) is a "marked region" that's present on
    one allele and absent on the other. Three primers tell the two alleles
    apart by PCR product size:
        F      : common forward primer  (upstream of the marked region, present on both alleles)
        R_in   : reverse primer picked *inside* the marked region -> only binds on the
                 allele that has that region
        R_out  : reverse primer picked *outside/downstream* of the marked region ->
                 binds on both alleles' underlying genomic sequence, but F+R_out only
                 falls in the amplifiable size range on the allele where the marked
                 region is absent (on the allele that has it, F+R_out spans the whole
                 region and is too large to amplify normally)

    Concretely, this same engine is used two ways:
        - Deletion/knockout (KO): the marked region is the deleted sequence, present
          in the WT allele and absent from the KO allele. Here R_in is called R_wt
          (WT-specific) and R_out is called R_ko (KO-specific) -- see
          design_3primer_genotyping() below.
        - HDR knock-in (KI): the marked region is the inserted cassette, present in
          the KI allele and absent from the WT allele -- the exact mirror image. Here
          R_in is called R_ki (KI-specific) and R_out is called R_wt (WT-specific) --
          see design_3primer_knockin_genotyping() further down, which reuses this same
          engine internally and relabels the output accordingly.

    Put all 3 primers in one tube and PCR (KO example):
        WT allele : F + R_wt  => product_wt  (short band)
        KO allele : F + R_ko  => product_ko  (shorter distance due to the deletion)
        (In the WT allele, F+R_ko spans the entire deletion length, so it's too
         large to amplify under normal PCR conditions)

    Constraints:
        product_min <= product_wt <= product_max   (default 200-1500bp)
        product_min <= product_ko <= product_max
        abs(product_wt - product_ko) >= min_size_diff (default 500bp, so the two
        bands are clearly distinguishable on a gel)

Genome-wide off-target (specificity) check (optional):
    primer3-py only checks self-dimer/hairpin/pair-dimer; it does NOT check whether
    a primer also binds elsewhere in the genome (that's what NCBI Primer-BLAST does
    on top of primer3). To do this locally you need:
        1) NCBI BLAST+ installed (Mac: `brew install blast`)
        2) genome fasta files ready (by default this looks in the genomes/ subfolder
           next to this script for
           GCF_000001405.40_GRCh38.p14_genomic.fna (hg38) and
           GCA_000001405.14_GRCh37.p13_genomic.fna (hg19))
    With both of these in place, design_3primer_genotyping() automatically aligns
    every candidate primer against the genome with blastn (-task blastn-short), and
    discards any candidate with 2 or more near-perfect hits (= suspected off-target).
    The BLAST DB is only built once (via makeblastdb); after that the cached DB is
    reused. Pass genome_fastas=[] to turn this check off entirely.

See the `if __name__ == "__main__":` block at the bottom of the file for usage examples.

Knock-in (HDR, large cassette insertion) genotyping:
    The functions above are all written for a deletion/knockout. For an HDR
    knock-in that inserts a large cassette (a reporter, tag, selection marker,
    etc.), use `expand_flanks_from_genome_using_arms()` (give it your left and
    right HDR homology arms; it locates both in the reference genome and
    figures out automatically whatever -- if anything -- sits between them),
    then `design_3primer_knockin_genotyping()` and `print_knockin_design()`
    near the bottom of this file. See `run_design_KI.py` for a
    ready-to-fill-in example, and `run_design_KO.py` (or `run_design.py`) for
    the deletion/knockout example.
"""

import os
import re
import shutil
import subprocess

import primer3

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Genome-wide specificity check (local BLAST, hg38/hg19, etc.)
# ---------------------------------------------------------------------------
#
# primer3-py only checks "does this primer form a bad structure with itself/its
# partner", not "does this primer also bind elsewhere in the genome" (off-target).
# That's normally what NCBI Primer-BLAST adds on top of primer3 output. To do it
# locally without internet access you need NCBI BLAST+ (blastn, makeblastdb) and
# a genome fasta.

DEFAULT_GENOME_FASTAS = [
    os.path.join(BASE_DIR, "genomes", "GCF_000001405.40_GRCh38.p14_genomic.fna"),  # hg38 / GRCh38
    os.path.join(BASE_DIR, "genomes", "GCA_000001405.14_GRCh37.p13_genomic.fna"),  # hg19 / GRCh37
]

# Passing a short name like genome_build="hg19" resolves to the actual file path
# through this mapping.
GENOME_BUILD_PATHS = {
    "hg38": os.path.join(BASE_DIR, "genomes", "GCF_000001405.40_GRCh38.p14_genomic.fna"),
    "grch38": os.path.join(BASE_DIR, "genomes", "GCF_000001405.40_GRCh38.p14_genomic.fna"),
    "hg19": os.path.join(BASE_DIR, "genomes", "GCA_000001405.14_GRCh37.p13_genomic.fna"),
    "grch37": os.path.join(BASE_DIR, "genomes", "GCA_000001405.14_GRCh37.p13_genomic.fna"),
}


def resolve_genome_build(genome_build):
    """'hg19' / 'hg38' / 'grch37' / 'grch38' / 'both' (case-insensitive) -> list of
    genome fasta paths.

    Raises:
        ValueError: unknown build name, or the corresponding genome fasta file
            doesn't exist
    """
    key = genome_build.strip().lower()
    if key == "both":
        paths = list(GENOME_BUILD_PATHS.values())
        # de-duplicate, keeping hg38 then hg19 order
        seen = []
        for p in [GENOME_BUILD_PATHS["hg38"], GENOME_BUILD_PATHS["hg19"]]:
            if p not in seen:
                seen.append(p)
        paths = seen
    elif key in GENOME_BUILD_PATHS:
        paths = [GENOME_BUILD_PATHS[key]]
    else:
        raise ValueError(
            f"genome_build={genome_build!r}: unrecognized value. "
            f"Use one of 'hg19', 'hg38', 'grch37', 'grch38', 'both'."
        )

    missing = [p for p in paths if not os.path.isfile(p)]
    if missing:
        raise FileNotFoundError(
            f"Could not find the fasta file(s) for genome_build={genome_build!r}: "
            f"{missing}\nCheck that the files exist in the genomes/ folder."
        )
    return paths


def _resolve_deletion_coords(full_seq, deletion_seq, del_start, del_end):
    """Determine 0-based (del_start, del_end) [python slice convention] from either
    deletion_seq (a sequence) or (del_start, del_end) (coordinates).

    Note: del_start == del_end (a zero-width point) is allowed and is used by
    the knock-in helpers below to represent a pure insertion point (nothing
    removed). A true deletion should have del_start < del_end."""
    if del_start is None or del_end is None:
        if deletion_seq is None:
            raise ValueError("You must supply either deletion_seq or (del_start, del_end)")
        deletion_seq = deletion_seq.strip()
        if deletion_seq == "":
            raise ValueError("deletion_seq is empty. For a pure insertion point with nothing "
                              "removed, pass del_start/del_end coordinates directly instead "
                              "(with del_start == del_end).")
        idx = full_seq.find(deletion_seq)
        if idx == -1:
            raise ValueError("deletion_seq was not found as an exact match inside full_seq. "
                              "Pass coordinates (del_start/del_end) directly instead.")
        del_start, del_end = idx, idx + len(deletion_seq)

    if not (0 <= del_start <= del_end <= len(full_seq)):
        raise ValueError(f"del_start/del_end coordinates are invalid: {del_start}, {del_end} "
                          f"(full_seq length={len(full_seq)})")
    return del_start, del_end


INSERTION_MARKER = "^"


def resolve_knockin_insertion(unedited_seq, replaced_seq):
    """Figure out where a cassette goes into an unedited (WT) sequence, for a
    knock-in. Returns (clean_unedited_seq, insert_start, insert_end) as 0-based
    python-slice coordinates marking the span in clean_unedited_seq that the
    cassette will replace.

    replaced_seq:
        - An exact substring of unedited_seq that gets removed/replaced by the
          cassette (e.g. a STOP codon being removed while inserting a tag).
          Its location is found automatically, same as deletion_seq for a
          knockout.
        - None, "", or "NA" (case-insensitive) for a *pure insertion* with
          nothing removed. In that case, mark the exact insertion point
          directly inside unedited_seq with a single '^' character, e.g.
          "...ATGAAG^CTGGAT..." — it's stripped out automatically and its
          position becomes the (zero-width) insertion point.
    """
    unedited_seq = unedited_seq.strip()

    is_pure_insertion = replaced_seq is None or replaced_seq.strip().upper() in ("", "NA")

    if not is_pure_insertion:
        replaced_seq = replaced_seq.strip()
        idx = unedited_seq.find(replaced_seq)
        if idx == -1:
            raise ValueError("replaced_seq was not found as an exact match inside "
                              "unedited_seq. Pass 'NA' instead if nothing is being "
                              "replaced (pure insertion).")
        return unedited_seq, idx, idx + len(replaced_seq)

    n_markers = unedited_seq.count(INSERTION_MARKER)
    if n_markers != 1:
        raise ValueError(
            f"replaced_seq is 'NA' (pure insertion), so mark the exact insertion "
            f"point inside unedited_seq with a single '{INSERTION_MARKER}' character "
            f"(e.g. '...ATGAAG{INSERTION_MARKER}CTGGAT...'). Found {n_markers} "
            f"marker(s) instead of exactly 1."
        )
    idx = unedited_seq.index(INSERTION_MARKER)
    clean_seq = unedited_seq[:idx] + unedited_seq[idx + 1:]
    return clean_seq, idx, idx


def _check_blast_binary(path, name):
    """Supports both a bare executable name on PATH (e.g. 'blastn') and a directly
    specified absolute path. Returns the resolved executable path string."""
    found = shutil.which(path)
    if found:
        return found
    if os.path.isfile(path) and os.access(path, os.X_OK):
        return path
    raise FileNotFoundError(
        f"Could not find the {name} executable: {path}\n"
        f"NCBI BLAST+ must be installed. Mac: `brew install blast`\n"
        f"(see https://ncbi.github.io/blast/)"
    )


def ensure_blast_db(genome_fasta, makeblastdb_bin="makeblastdb"):
    """Build a BLAST nucleotide DB for genome_fasta if one doesn't already exist
    (one time only; can take minutes to tens of minutes depending on genome size).
    If it already exists, it's reused as-is."""
    if not os.path.isfile(genome_fasta):
        raise FileNotFoundError(f"Could not find the genome fasta file: {genome_fasta}")
    makeblastdb_bin = _check_blast_binary(makeblastdb_bin, "makeblastdb")
    db_markers = [genome_fasta + ext for ext in (".nin", ".ndb", ".00.nin")]
    if not any(os.path.exists(m) for m in db_markers):
        print(f"[genome specificity] No BLAST DB found for '{os.path.basename(genome_fasta)}', "
              f"building one now (one time only; may take a while depending on genome size)...")
        subprocess.run([makeblastdb_bin, "-in", genome_fasta, "-dbtype", "nucl",
                         "-parse_seqids"], check=True)
    return genome_fasta


def _run_blastn_short(seq, db, blastn_bin="blastn", evalue=1000, word_size=7,
                       perc_identity=80):
    """Align a single short primer (18-30nt) against the genome DB using blastn-short."""
    blastn_bin = _check_blast_binary(blastn_bin, "blastn")
    query = f">primer\n{seq}\n"
    cmd = [blastn_bin, "-task", "blastn-short", "-db", db,
           "-word_size", str(word_size), "-evalue", str(evalue),
           "-perc_identity", str(perc_identity),
           "-outfmt", "6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"]
    proc = subprocess.run(cmd, input=query, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"blastn failed:\n{proc.stderr}")
    hits = []
    for line in proc.stdout.strip().splitlines():
        if not line:
            continue
        p = line.split("\t")
        hits.append(dict(sseqid=p[0], pident=float(p[1]), length=int(p[2]),
                          mismatch=int(p[3]), gapopen=int(p[4]),
                          qstart=int(p[5]), qend=int(p[6]),
                          sstart=int(p[7]), send=int(p[8]),
                          evalue=float(p[9]), bitscore=float(p[10])))
    return hits


_SPECIFICITY_CACHE = {}


def check_genome_specificity(primer_seq, genome_fastas,
                              blastn_bin="blastn", makeblastdb_bin="makeblastdb",
                              min_pident=90.0, min_coverage=0.9, require_3prime=True,
                              evalue=1000, word_size=7):
    """Check how many places in each of genome_fastas primer_seq binds "almost
    perfectly".

    A single hit (= the intended target site) is considered specific; 2 or more
    hits indicates off-target risk. (This function doesn't know which hit is the
    "intended" one, so n_hits==1 is usually assumed safe, while n_hits>=2 needs a
    human to check which hit is the real target.)

    Returns: {genome_fasta_path: {"n_hits": int, "hits": [...]}}
    """
    if isinstance(genome_fastas, str):
        genome_fastas = [genome_fastas]

    primer_seq = primer_seq.upper()
    result = {}
    for genome_fasta in genome_fastas:
        cache_key = (primer_seq, genome_fasta, min_pident, min_coverage, require_3prime)
        if cache_key in _SPECIFICITY_CACHE:
            result[genome_fasta] = _SPECIFICITY_CACHE[cache_key]
            continue

        db = ensure_blast_db(genome_fasta, makeblastdb_bin=makeblastdb_bin)
        hits = _run_blastn_short(primer_seq, db, blastn_bin=blastn_bin,
                                  evalue=evalue, word_size=word_size)

        good_hits = []
        for h in hits:
            if h["pident"] < min_pident:
                continue
            if h["length"] < min_coverage * len(primer_seq):
                continue
            # The 3' end (where extension actually happens) must be included in the
            # alignment for this to count as a genuinely risky hit.
            if require_3prime and h["qend"] < len(primer_seq) - 1:
                continue
            good_hits.append(h)

        entry = dict(n_hits=len(good_hits), hits=good_hits)
        _SPECIFICITY_CACHE[cache_key] = entry
        result[genome_fasta] = entry

    return result


def is_primer_specific(primer_seq, genome_fastas, max_hits=1, **kwargs):
    """True if the number of near-perfect hits is <= max_hits in every genome_fastas."""
    spec = check_genome_specificity(primer_seq, genome_fastas, **kwargs)
    return all(v["n_hits"] <= max_hits for v in spec.values())


# ---------------------------------------------------------------------------
# Automatically fetch generous flanking sequence around the deletion from the genome
# ---------------------------------------------------------------------------
#
# The full_seq a user supplies usually only has sequence immediately around the
# deletion, which can make the search space for the KO-specific reverse primer
# (located "downstream" of where the deletion ends) too tight. This feature
# (1) BLASTs full_seq against the genome to find its exact location, then
# (2) attaches flank_kb (default 50kb) of sequence on each side of the deletion,
# producing a new full_seq/del_start/del_end. blastdbcmd (bundled with BLAST+) is
# used to cut the region directly out of the genome fasta.

def _run_blastn_full(seq, db, blastn_bin="blastn", task="megablast", evalue=1e-20):
    """BLAST the entire full_seq (hundreds to thousands of bp) against the genome DB
    to locate it."""
    blastn_bin = _check_blast_binary(blastn_bin, "blastn")
    query = f">query\n{seq}\n"
    cmd = [blastn_bin, "-task", task, "-db", db, "-evalue", str(evalue),
           "-outfmt", "6 sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen"]
    proc = subprocess.run(cmd, input=query, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"blastn failed:\n{proc.stderr}")
    hits = []
    for line in proc.stdout.strip().splitlines():
        if not line:
            continue
        p = line.split("\t")
        hits.append(dict(sseqid=p[0], pident=float(p[1]), length=int(p[2]),
                          mismatch=int(p[3]), gapopen=int(p[4]),
                          qstart=int(p[5]), qend=int(p[6]),
                          sstart=int(p[7]), send=int(p[8]),
                          evalue=float(p[9]), bitscore=float(p[10]),
                          slen=int(p[11])))
    return hits


def _extract_genome_range(db, sseqid, start_1based, end_1based, blastdbcmd_bin="blastdbcmd"):
    """Extract the [start_1based, end_1based] range (1-based, inclusive, always
    plus/forward strand) directly from the genome fasta (BLAST DB)."""
    blastdbcmd_bin = _check_blast_binary(blastdbcmd_bin, "blastdbcmd")
    cmd = [blastdbcmd_bin, "-db", db, "-entry", sseqid,
           "-range", f"{start_1based}-{end_1based}", "-strand", "plus"]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"blastdbcmd failed:\n{proc.stderr}")
    seq_lines = [ln for ln in proc.stdout.strip().splitlines() if not ln.startswith(">")]
    return "".join(seq_lines).upper()


def _locate_full_seq_genomic_span(seq, db, blastn_bin, min_pident, min_coverage, label):
    """BLAST a whole sequence (e.g. one HDR homology arm) against the genome DB
    and return the ascending plus-strand coordinate range it occupies, along
    with which strand it aligned to.

    Returns: dict(sseqid, strand_plus, genome_start, genome_end, slen), where
    genome_start <= genome_end always (the plus-strand coordinate interval the
    sequence occupies), regardless of which strand it aligns to.
    """
    hits = _run_blastn_full(seq, db, blastn_bin=blastn_bin)
    candidates = []
    for h in hits:
        if h["pident"] < min_pident:
            continue
        if h["length"] < min_coverage * len(seq):
            continue
        if h["gapopen"] != 0:
            continue
        candidates.append(h)
    if not candidates:
        raise RuntimeError(
            f"Could not reliably locate {label} in the genome (insufficient "
            f"identity/coverage, or a gapped alignment). Check that genome_build "
            f"matches where this sequence came from."
        )
    best = max(candidates, key=lambda h: h["bitscore"])
    strand_plus = best["sstart"] < best["send"]
    if strand_plus:
        genome_start = best["sstart"] - (best["qstart"] - 1)
        genome_end = best["sstart"] + (len(seq) - 1 - (best["qstart"] - 1))
    else:
        genome_start = best["sstart"] - (len(seq) - 1 - (best["qstart"] - 1))
        genome_end = best["sstart"] + (best["qstart"] - 1)
    return dict(sseqid=best["sseqid"], strand_plus=strand_plus,
                genome_start=genome_start, genome_end=genome_end, slen=best["slen"])


def expand_flanks_from_genome_using_arms(
    left_arm,
    right_arm,
    genome_build=None,
    genome_fastas=None,
    flank_kb=50,
    min_pident=98.0,
    min_coverage=0.95,
    blastn_bin="blastn",
    makeblastdb_bin="makeblastdb",
    blastdbcmd_bin="blastdbcmd",
):
    """Locate the left and right HDR homology arms in the reference genome, and
    fetch a single contiguous WT-locus sequence spanning from flank_kb
    upstream of the left arm to flank_kb downstream of the right arm.

    This is the simplest way to set up a knock-in design: whatever sequence
    naturally sits between the two arms in the genome (nothing, for a pure
    insertion with the arms directly adjacent -- or a real stretch, for a
    design that replaces/removes a short region) is figured out automatically
    from the genome itself. You don't need to separately specify a
    "replaced_seq" or mark an insertion point by hand.

    Both arms must be given in plus-strand (reference/forward) orientation;
    if your HDR design is on the minus strand, reverse-complement both arms
    first (this function currently only supports plus-strand loci).

    Returns:
        (wt_locus_seq, insert_start, insert_end, info)
        wt_locus_seq   : the fetched genomic sequence (plus strand)
        insert_start, insert_end : 0-based python-slice coordinates in
            wt_locus_seq marking the span between the two arms, where the
            cassette goes in (insert_start == insert_end for a pure insertion
            with the arms directly adjacent in the genome)
        info : chrom, each arm's genomic coordinates, the gap length between
            them, and the upstream/downstream flank actually obtained
    """
    left_arm = left_arm.strip().upper()
    right_arm = right_arm.strip().upper()

    if genome_fastas is None:
        if genome_build is None:
            raise ValueError(
                "expand_flanks_from_genome_using_arms() needs exactly one genome: "
                "specify genome_build='hg19' or 'hg38' (or 'grch37'/'grch38')."
            )
        genome_fastas = resolve_genome_build(genome_build)
    if isinstance(genome_fastas, str):
        genome_fastas = [genome_fastas]
    if len(genome_fastas) != 1:
        raise ValueError(
            f"expand_flanks_from_genome_using_arms() requires exactly one genome "
            f"(received {len(genome_fastas)}). Use a single value like genome_build='hg19'."
        )
    genome_fasta = genome_fastas[0]
    db = ensure_blast_db(genome_fasta, makeblastdb_bin=makeblastdb_bin)

    print(f"[genome context] Locating the left and right homology arms in "
          f"'{os.path.basename(genome_fasta)}' (blastn)...")
    left_loc = _locate_full_seq_genomic_span(left_arm, db, blastn_bin, min_pident,
                                              min_coverage, "left_arm")
    right_loc = _locate_full_seq_genomic_span(right_arm, db, blastn_bin, min_pident,
                                               min_coverage, "right_arm")

    if left_loc["sseqid"] != right_loc["sseqid"]:
        raise RuntimeError(
            f"left_arm and right_arm aligned to different chromosomes/contigs "
            f"({left_loc['sseqid']} vs {right_loc['sseqid']}). They should come "
            f"from the same locus."
        )
    if (not left_loc["strand_plus"]) or (not right_loc["strand_plus"]):
        raise RuntimeError(
            "One or both arms aligned to the minus strand, which isn't currently "
            "supported. Reverse-complement both left_arm and right_arm (so they "
            "align to the plus strand of the reference) and try again."
        )
    if left_loc["genome_end"] >= right_loc["genome_start"]:
        raise RuntimeError(
            f"left_arm (genomic {left_loc['genome_start']}-{left_loc['genome_end']}) "
            f"does not fall entirely upstream of right_arm (genomic "
            f"{right_loc['genome_start']}-{right_loc['genome_end']}). Check that the "
            f"arms are in the correct left/right order and don't overlap."
        )

    slen = left_loc["slen"]
    ext_start = max(1, left_loc["genome_start"] - flank_kb * 1000)
    ext_end = min(slen, right_loc["genome_end"] + flank_kb * 1000)
    gap_len = right_loc["genome_start"] - left_loc["genome_end"] - 1

    print(f"[genome context] Located: {left_loc['sseqid']} "
          f"left_arm={left_loc['genome_start']}-{left_loc['genome_end']}, "
          f"right_arm={right_loc['genome_start']}-{right_loc['genome_end']} "
          f"(gap between arms: {gap_len:,}bp) -> extracting {ext_start}-{ext_end} "
          f"({ext_end - ext_start + 1:,}bp)")

    wt_locus_seq = _extract_genome_range(db, left_loc["sseqid"], ext_start, ext_end,
                                          blastdbcmd_bin=blastdbcmd_bin)

    insert_start = left_loc["genome_end"] - ext_start + 1
    insert_end = right_loc["genome_start"] - ext_start

    # Defensive sanity check: the arms should reappear at exactly these
    # positions in what we just extracted. If not, something is off with the
    # coordinate math or the underlying BLAST alignment.
    if wt_locus_seq[insert_start - len(left_arm):insert_start] != left_arm:
        raise RuntimeError("Internal consistency check failed: the extracted sequence "
                            "doesn't reproduce left_arm at the expected position. "
                            "This may indicate mismatches near the arm boundary.")
    if wt_locus_seq[insert_end:insert_end + len(right_arm)] != right_arm:
        raise RuntimeError("Internal consistency check failed: the extracted sequence "
                            "doesn't reproduce right_arm at the expected position. "
                            "This may indicate mismatches near the arm boundary.")

    info = dict(
        chrom=left_loc["sseqid"], strand="plus",
        left_arm_genomic=(left_loc["genome_start"], left_loc["genome_end"]),
        right_arm_genomic=(right_loc["genome_start"], right_loc["genome_end"]),
        gap_len=gap_len,
        ext_start=ext_start, ext_end=ext_end,
        upstream_flank_bp=insert_start - len(left_arm),
        downstream_flank_bp=len(wt_locus_seq) - insert_end - len(right_arm),
    )
    return wt_locus_seq, insert_start, insert_end, info


def expand_flanks_from_genome(
    full_seq,
    deletion_seq=None,
    del_start=None,
    del_end=None,
    genome_build=None,
    genome_fastas=None,
    flank_kb=50,
    min_pident=98.0,
    min_coverage=0.95,
    blastn_bin="blastn",
    makeblastdb_bin="makeblastdb",
    blastdbcmd_bin="blastdbcmd",
):
    """Fetch a new sequence with flank_kb (default 50kb) attached on each side of
    the deletion, taken from the genome.

    How it works:
        1) Align full_seq against the genome (BLAST DB) to find the exact genomic
           coordinates (identity >= min_pident %, coverage >= min_coverage, only
           gapless alignments accepted, and both deletion boundaries must fall
           within the aligned range)
        2) Compute the genomic coordinates of the deletion boundaries
        3) Use blastdbcmd to cut out the whole
           (deletion_start - flank_kb*1000) to (deletion_end + flank_kb*1000) range
           as the new full_seq (automatically clamped if it runs past the end of
           the chromosome)

    Args:
        genome_build: must resolve to exactly one genome ("hg19"/"hg38"/"grch37"/
                      "grch38"). "both" cannot be used here since we need to know
                      exactly which genome to cut from.
        genome_fastas: use this instead of genome_build to specify a single file
                       path directly.

    Returns:
        A (new_full_seq, new_del_start, new_del_end, info) tuple.
        info contains the mapped genomic coordinates, the chromosome accession,
        the strand, and the upstream/downstream flank lengths actually obtained
        (may be smaller than the requested flank_kb if clamped near the end of
        the chromosome).
    """
    full_seq = full_seq.strip()
    del_start, del_end = _resolve_deletion_coords(full_seq, deletion_seq, del_start, del_end)

    if genome_fastas is None:
        if genome_build is None:
            raise ValueError(
                "expand_flanks_from_genome() needs exactly one genome to pull from: "
                "specify genome_build='hg19' or 'hg38' (or 'grch37'/'grch38')."
            )
        genome_fastas = resolve_genome_build(genome_build)
    if isinstance(genome_fastas, str):
        genome_fastas = [genome_fastas]
    if len(genome_fastas) != 1:
        raise ValueError(
            f"expand_flanks_from_genome() requires exactly one genome "
            f"(received {len(genome_fastas)}). Use a single value like genome_build='hg19'."
        )
    genome_fasta = genome_fastas[0]
    db = ensure_blast_db(genome_fasta, makeblastdb_bin=makeblastdb_bin)

    print(f"[genome context] Locating this sequence in '{os.path.basename(genome_fasta)}' "
          f"(blastn)...")
    hits = _run_blastn_full(full_seq, db, blastn_bin=blastn_bin)

    candidates = []
    for h in hits:
        if h["pident"] < min_pident:
            continue
        if h["length"] < min_coverage * len(full_seq):
            continue
        if h["gapopen"] != 0:
            continue
        # Both deletion boundaries (1-based: del_start+1, del_end) must fall within
        # the aligned range for the coordinate conversion to be accurate.
        if not (h["qstart"] <= del_start + 1 <= h["qend"]):
            continue
        if not (h["qstart"] <= del_end <= h["qend"]):
            continue
        candidates.append(h)

    if not candidates:
        raise RuntimeError(
            "Could not reliably locate this sequence in the genome "
            "(insufficient identity/coverage, a gapped alignment, or the deletion "
            "boundaries fall outside the aligned range). Check whether genome_build "
            "actually matches where this sequence came from (hg19 vs hg38)."
        )

    best = max(candidates, key=lambda h: h["bitscore"])
    strand_plus = best["sstart"] < best["send"]

    if strand_plus:
        genome_del_start = best["sstart"] + (del_start - (best["qstart"] - 1))
        genome_del_end = best["sstart"] + (del_end - 1 - (best["qstart"] - 1))
    else:
        genome_del_start = best["sstart"] - (del_end - 1 - (best["qstart"] - 1))
        genome_del_end = best["sstart"] - (del_start - (best["qstart"] - 1))

    slen = best["slen"]
    ext_start = max(1, genome_del_start - flank_kb * 1000)
    ext_end = min(slen, genome_del_end + flank_kb * 1000)

    print(f"[genome context] Located at: {best['sseqid']}:{genome_del_start}-{genome_del_end} "
          f"({'plus' if strand_plus else 'minus'} strand, identity={best['pident']:.1f}%) -> "
          f"extracting {ext_start}-{ext_end} ({ext_end - ext_start + 1:,}bp)")

    new_full_seq = _extract_genome_range(db, best["sseqid"], ext_start, ext_end,
                                          blastdbcmd_bin=blastdbcmd_bin)
    new_del_start = genome_del_start - ext_start
    new_del_end = genome_del_end - ext_start + 1

    info = dict(
        chrom=best["sseqid"], strand="plus" if strand_plus else "minus",
        genome_del_start=genome_del_start, genome_del_end=genome_del_end,
        ext_start=ext_start, ext_end=ext_end,
        upstream_flank_bp=new_del_start,
        downstream_flank_bp=len(new_full_seq) - new_del_end,
    )
    return new_full_seq, new_del_start, new_del_end, info


# ---------------------------------------------------------------------------
# Get candidate primers via primer3-py (PRIMER_TASK=pick_primer_list)
# ---------------------------------------------------------------------------

def _parse_coord(val, is_left):
    """Convert a PRIMER_LEFT_i / PRIMER_RIGHT_i coordinate value into (start, end)
    [python slice convention]. Depending on the primer3-py version, val may be the
    string '1432,20' or the tuple (1432, 20)."""
    if isinstance(val, (tuple, list)):
        pos, length = int(val[0]), int(val[1])
    elif isinstance(val, str):
        a, b = val.split(",")
        pos, length = int(a), int(b)
    else:
        raise TypeError(f"Unexpected coordinate value: {val!r}")

    if is_left:
        return pos, pos + length, length
    else:
        end = pos + 1
        start = pos - length + 1
        return start, end, length


def get_candidates(side, full_seq, region_start, region_len,
                    len_range=(18, 28), tm_range=(59, 63), tm_opt=60.0,
                    gc_range=(40, 60), gc_opt=50.0,
                    max_self_any_th=47.0, max_self_end_th=47.0, max_hairpin_th=47.0,
                    num_return=300, lowercase_masking=True):
    """side: 'left' (forward candidates) or 'right' (reverse candidates).
    region_start, region_len: the 0-based window within full_seq where this primer
    is allowed to sit (SEQUENCE_INCLUDED_REGION)."""
    is_left = (side == "left")

    seq_args = {
        "SEQUENCE_ID": "genotyping_seq",
        "SEQUENCE_TEMPLATE": full_seq,
        "SEQUENCE_INCLUDED_REGION": (region_start, region_len),
    }
    global_args = {
        "PRIMER_TASK": "pick_primer_list",
        "PRIMER_PICK_LEFT_PRIMER": 1 if is_left else 0,
        "PRIMER_PICK_RIGHT_PRIMER": 0 if is_left else 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_MIN_SIZE": len_range[0],
        "PRIMER_MAX_SIZE": len_range[1],
        "PRIMER_OPT_SIZE": sum(len_range) // 2,
        "PRIMER_MIN_TM": tm_range[0],
        "PRIMER_MAX_TM": tm_range[1],
        "PRIMER_OPT_TM": tm_opt,
        "PRIMER_MIN_GC": gc_range[0],
        "PRIMER_MAX_GC": gc_range[1],
        "PRIMER_OPT_GC_PERCENT": gc_opt,
        "PRIMER_MAX_SELF_ANY_TH": max_self_any_th,
        "PRIMER_MAX_SELF_END_TH": max_self_end_th,
        "PRIMER_MAX_HAIRPIN_TH": max_hairpin_th,
        "PRIMER_LOWERCASE_MASKING": 1 if lowercase_masking else 0,
    }

    out = primer3.bindings.design_primers(seq_args, global_args)
    if "PRIMER_ERROR" in out:
        raise RuntimeError(f"primer3 error: {out['PRIMER_ERROR']}")

    prefix = "PRIMER_LEFT_" if is_left else "PRIMER_RIGHT_"
    n_key = "PRIMER_LEFT_NUM_RETURNED" if is_left else "PRIMER_RIGHT_NUM_RETURNED"
    n = int(out.get(n_key, 0))

    cands = []
    for i in range(n):
        seq = str(out[f"{prefix}{i}_SEQUENCE"]).upper()
        start, end, length = _parse_coord(out[f"{prefix}{i}"], is_left)
        cands.append(dict(
            seq=seq, start=start, end=end, len=length,
            tm=float(out[f"{prefix}{i}_TM"]),
            gc=float(out[f"{prefix}{i}_GC_PERCENT"]),
            self_any=float(out.get(f"{prefix}{i}_SELF_ANY_TH", 0.0)),
            self_end=float(out.get(f"{prefix}{i}_SELF_END_TH", 0.0)),
            hairpin=float(out.get(f"{prefix}{i}_HAIRPIN_TH", 0.0)),
        ))
    return cands


def pair_compl_any(seq1, seq2):
    """Overall complementarity between two primers (corresponds to the Tm scale of
    PRIMER_PAIR_COMPL_ANY_TH)."""
    return primer3.bindings.calc_heterodimer(seq1, seq2).tm


def pair_compl_end(seq1, seq2):
    """3' end complementarity between two primers (corresponds to
    PRIMER_PAIR_COMPL_END_TH). Returns the riskier of the two directions (seq1's
    3' end vs seq2, and seq2's 3' end vs seq1)."""
    a = primer3.bindings.calc_end_stability(seq1, seq2).tm
    b = primer3.bindings.calc_end_stability(seq2, seq1).tm
    return max(a, b)


# ---------------------------------------------------------------------------
# Main design function
# ---------------------------------------------------------------------------

def design_3primer_genotyping(
    full_seq,
    deletion_seq=None,
    del_start=None,
    del_end=None,
    product_min=200,
    product_max=1500,
    min_size_diff=500,
    primer_len_range=(18, 28),
    tm_range=(59, 63),
    tm_opt=60.0,
    gc_range=(40, 60),
    gc_opt=50.0,
    max_search_upstream=1600,
    max_search_downstream=1600,
    max_search_into_deletion=1600,
    max_self_any_th=47.0,
    max_self_end_th=47.0,
    max_hairpin_th=47.0,
    max_pair_compl_any_th=47.0,
    max_pair_compl_end_th=47.0,
    max_tm_spread=3.0,
    top_n=5,
    candidate_pool_size=150,
    # --- genome-wide off-target specificity (local BLAST) ---
    genome_build=None,                # "hg19" / "hg38" / "grch37" / "grch38" / "both"
                                       # Supplying this auto-fills genome_fastas (simplest option)
    genome_fastas=None,               # None (with genome_build also None) => auto-discover both in genomes/
                                       # []   => turn off specificity check entirely
                                       # ["/path/to/genome.fa", ...] => specify file path(s) directly
    max_genome_hits=1,
    specificity_min_pident=90.0,
    specificity_min_coverage=0.9,
    blastn_bin="blastn",
    makeblastdb_bin="makeblastdb",
):
    """Design a 3-primer genotyping set. Returns a list of dicts (up to top_n),
    sorted ascending by score (lower = better). Each dict holds the F/R_wt/R_ko
    primer info, product sizes, pair-dimer scores, etc."""
    full_seq = full_seq.strip()

    # ---- 1) Determine deletion coordinates ----
    del_start, del_end = _resolve_deletion_coords(full_seq, deletion_seq, del_start, del_end)

    common_kw = dict(len_range=primer_len_range, tm_range=tm_range, tm_opt=tm_opt,
                      gc_range=gc_range, gc_opt=gc_opt,
                      max_self_any_th=max_self_any_th, max_self_end_th=max_self_end_th,
                      max_hairpin_th=max_hairpin_th, num_return=candidate_pool_size)

    # ---- 2) Get candidate primer pools for the three regions (real primer3-py calls) ----
    f_region_start = max(0, del_start - max_search_upstream)
    f_region_len = del_start - f_region_start
    fwd_cands = get_candidates("left", full_seq, f_region_start, f_region_len, **common_kw)

    rwt_region_start = del_start
    rwt_region_len = min(max_search_into_deletion, del_end - del_start)
    rwt_cands = get_candidates("right", full_seq, rwt_region_start, rwt_region_len, **common_kw)

    rko_region_start = del_end
    rko_region_len = min(max_search_downstream, len(full_seq) - del_end)
    rko_cands = get_candidates("right", full_seq, rko_region_start, rko_region_len, **common_kw)

    if not fwd_cands:
        raise RuntimeError("Could not find any common forward primer candidates. "
                            "Try increasing max_search_upstream or widening the tm/gc ranges.")
    if not rwt_cands:
        raise RuntimeError("Could not find any WT-specific reverse primer candidates. "
                            "The deletion region may be too short, or the constraints too tight.")
    if not rko_cands:
        raise RuntimeError("Could not find any KO-specific reverse primer candidates. "
                            "There may not be enough sequence past the end of the deletion "
                            "(downstream flank needed), or the constraints may be too tight.")

    pool_stats = dict(
        fwd_raw=len(fwd_cands), rwt_raw=len(rwt_cands), rko_raw=len(rko_cands),
    )

    def _topk(cands, k=40):
        return sorted(cands, key=lambda c: abs(c["tm"] - tm_opt))[:k]

    fwd_top = _topk(fwd_cands)
    rwt_top = _topk(rwt_cands)
    rko_top = _topk(rko_cands)
    pool_stats.update(
        fwd_topk=len(fwd_top), rwt_topk=len(rwt_top), rko_topk=len(rko_top),
    )

    # ---- 2b) Genome-wide off-target specificity check (local BLAST) ----
    if genome_fastas is None:
        if genome_build is not None:
            genome_fastas = resolve_genome_build(genome_build)
            print(f"[genome specificity] genome_build={genome_build!r} -> "
                  f"{[os.path.basename(f) for f in genome_fastas]}")
        else:
            genome_fastas = [f for f in DEFAULT_GENOME_FASTAS if os.path.isfile(f)]
            if genome_fastas:
                print(f"[genome specificity] Auto-discovered genome fasta: "
                      f"{[os.path.basename(f) for f in genome_fastas]}")
            else:
                print("[genome specificity] No genome fasta found, skipping the off-target "
                      "check. Pass genome_build='hg19'/'hg38' or the genome_fastas argument "
                      "directly to enable it.")
    elif not genome_fastas:
        print("[genome specificity] genome_fastas=[] -- genome-wide off-target checking is "
              "disabled for this run.")

    spec_kw = dict(blastn_bin=blastn_bin, makeblastdb_bin=makeblastdb_bin,
                    min_pident=specificity_min_pident, min_coverage=specificity_min_coverage)

    def _filter_specific(cands):
        if not genome_fastas:
            for c in cands:
                c["specificity"] = None
            return cands
        kept = []
        for c in cands:
            spec = check_genome_specificity(c["seq"], genome_fastas, **spec_kw)
            c["specificity"] = spec
            if all(v["n_hits"] <= max_genome_hits for v in spec.values()):
                kept.append(c)
        return kept

    fwd_top = _filter_specific(fwd_top)
    rwt_top = _filter_specific(rwt_top)
    rko_top = _filter_specific(rko_top)
    pool_stats.update(
        fwd_after_spec=len(fwd_top), rwt_after_spec=len(rwt_top), rko_after_spec=len(rko_top),
    )

    if genome_fastas and not (fwd_top and rwt_top and rko_top):
        print_diagnostics(pool_stats, None)
        raise RuntimeError(
            "No candidates remain after genome specificity filtering (too many "
            "off-targets, or the constraints are too tight). Try increasing "
            "max_genome_hits or candidate_pool_size."
        )

    # ---- 3) Evaluate combinations ----
    dimer_cache = {}

    def _pair_score(seq1, seq2, mode):
        key = (seq1, seq2, mode)
        if key not in dimer_cache:
            if mode == "ANY":
                dimer_cache[key] = pair_compl_any(seq1, seq2)
            else:
                dimer_cache[key] = pair_compl_end(seq1, seq2)
        return dimer_cache[key]

    combo_stats = dict(
        total_possible=len(fwd_top) * len(rwt_top) * len(rko_top),
        fail_product_wt=0, fail_product_ko=0, fail_size_diff=0,
        fail_tm_spread=0, fail_pair_any=0, fail_pair_end=0, passed=0,
    )

    results = []
    for f in fwd_top:
        for rwt in rwt_top:
            product_wt = rwt["end"] - f["start"]
            if not (product_min <= product_wt <= product_max):
                combo_stats["fail_product_wt"] += len(rko_top)
                continue
            for rko in rko_top:
                product_ko = (del_start - f["start"]) + (rko["end"] - del_end)
                if not (product_min <= product_ko <= product_max):
                    combo_stats["fail_product_ko"] += 1
                    continue
                size_diff = abs(product_wt - product_ko)
                if size_diff < min_size_diff:
                    combo_stats["fail_size_diff"] += 1
                    continue

                tms = [f["tm"], rwt["tm"], rko["tm"]]
                tm_spread = max(tms) - min(tms)
                if tm_spread > max_tm_spread:
                    combo_stats["fail_tm_spread"] += 1
                    continue

                pd_f_rwt = _pair_score(f["seq"], rwt["seq"], "ANY")
                if pd_f_rwt > max_pair_compl_any_th:
                    combo_stats["fail_pair_any"] += 1
                    continue
                pd_f_rko = _pair_score(f["seq"], rko["seq"], "ANY")
                if pd_f_rko > max_pair_compl_any_th:
                    combo_stats["fail_pair_any"] += 1
                    continue
                pd_rwt_rko = _pair_score(rwt["seq"], rko["seq"], "ANY")
                if pd_rwt_rko > max_pair_compl_any_th:
                    combo_stats["fail_pair_any"] += 1
                    continue

                pd_f_rwt_end = _pair_score(f["seq"], rwt["seq"], "END")
                pd_f_rko_end = _pair_score(f["seq"], rko["seq"], "END")
                pd_rwt_rko_end = _pair_score(rwt["seq"], rko["seq"], "END")
                if max(pd_f_rwt_end, pd_f_rko_end, pd_rwt_rko_end) > max_pair_compl_end_th:
                    combo_stats["fail_pair_end"] += 1
                    continue

                combo_stats["passed"] += 1
                score = (
                    tm_spread * 3.0
                    + (f["self_any"] + rwt["self_any"] + rko["self_any"]) * 0.02
                    + (f["self_end"] + rwt["self_end"] + rko["self_end"]) * 0.05
                    + (pd_f_rwt + pd_f_rko + pd_rwt_rko) * 0.05
                    + (pd_f_rwt_end + pd_f_rko_end + pd_rwt_rko_end) * 0.1
                    - size_diff * 0.005
                )

                results.append(dict(
                    F=f, R_wt=rwt, R_ko=rko,
                    product_wt=product_wt, product_ko=product_ko,
                    size_diff=size_diff, tm_spread=tm_spread,
                    pair_compl_any=dict(F_Rwt=pd_f_rwt, F_Rko=pd_f_rko, Rwt_Rko=pd_rwt_rko),
                    pair_compl_end=dict(F_Rwt=pd_f_rwt_end, F_Rko=pd_f_rko_end, Rwt_Rko=pd_rwt_rko_end),
                    del_start=del_start, del_end=del_end,
                    score=score,
                ))

    results.sort(key=lambda r: r["score"])

    if not results:
        print_diagnostics(pool_stats, combo_stats)

    return results[:top_n]


# ---------------------------------------------------------------------------
# Diagnostics (shows where/how many candidates were dropped when design fails)
# ---------------------------------------------------------------------------

def print_diagnostics(pool_stats, combo_stats):
    print("\n=== Diagnostics: why no design was found ===")
    print(f"[Candidate pools]")
    print(f"  Forward : raw={pool_stats['fwd_raw']:>5}  -> top-k={pool_stats['fwd_topk']:>4}"
          + (f"  -> passed genome-specificity={pool_stats['fwd_after_spec']:>4}"
             if 'fwd_after_spec' in pool_stats else "  (specificity not checked)"))
    print(f"  R_wt    : raw={pool_stats['rwt_raw']:>5}  -> top-k={pool_stats['rwt_topk']:>4}"
          + (f"  -> passed genome-specificity={pool_stats['rwt_after_spec']:>4}"
             if 'rwt_after_spec' in pool_stats else "  (specificity not checked)"))
    print(f"  R_ko    : raw={pool_stats['rko_raw']:>5}  -> top-k={pool_stats['rko_topk']:>4}"
          + (f"  -> passed genome-specificity={pool_stats['rko_after_spec']:>4}"
             if 'rko_after_spec' in pool_stats else "  (specificity not checked)"))

    if combo_stats is None:
        print("\n(One of the candidate pools was emptied out during genome specificity "
              "filtering, so the combinatorial search never even started.)")
        return

    total = combo_stats["total_possible"]
    print(f"\n[Combinatorial search] Combinations attempted: {total:,}")
    if total == 0:
        print("  (One of the candidate pools is empty, so 0 combinations were tried.)")
        return

    def _pct(n):
        return f"{n:>6,} ({100*n/total:5.1f}%)" if total else f"{n:>6,}"

    print(f"  product_wt out of range (e.g. 200-1500)  : {_pct(combo_stats['fail_product_wt'])}")
    print(f"  product_ko out of range                  : {_pct(combo_stats['fail_product_ko'])}")
    print(f"  size_diff < min_size_diff                 : {_pct(combo_stats['fail_size_diff'])}")
    print(f"  Tm spread > max_tm_spread                 : {_pct(combo_stats['fail_tm_spread'])}")
    print(f"  pair complementarity (ANY) exceeded        : {_pct(combo_stats['fail_pair_any'])}")
    print(f"  pair complementarity (3' END) exceeded     : {_pct(combo_stats['fail_pair_end'])}")
    print(f"  passed                                    : {_pct(combo_stats['passed'])}")
    print()

    # Point out whichever stage filtered out the most combinations
    reasons = {
        "product_wt range": combo_stats["fail_product_wt"],
        "product_ko range": combo_stats["fail_product_ko"],
        "size_diff too small": combo_stats["fail_size_diff"],
        "Tm spread exceeded": combo_stats["fail_tm_spread"],
        "pair ANY exceeded": combo_stats["fail_pair_any"],
        "pair END exceeded": combo_stats["fail_pair_end"],
    }
    top_reason = max(reasons, key=reasons.get)
    if reasons[top_reason] > 0:
        print(f"-> Biggest bottleneck: the '{top_reason}' stage. Consider relaxing this "
              f"constraint, or growing the candidate pool itself by widening the flanks "
              f"with expand_flanks_from_genome().\n")


# ---------------------------------------------------------------------------
# Result printing
# ---------------------------------------------------------------------------

def _spec_str(cand):
    spec = cand.get("specificity")
    if not spec:
        return "specificity=not checked"
    parts = []
    for genome_path, v in spec.items():
        parts.append(f"{os.path.basename(genome_path)}:{v['n_hits']}hit(s)")
    return "genome_hits[" + ", ".join(parts) + "]"


def print_design(d, rank=None):
    print(f"=== Design{' #' + str(rank) if rank else ''} (score={d['score']:.2f}) ===")
    print(f"Deletion coordinates (0-based, python slice): {d['del_start']}-{d['del_end']} "
          f"(length {d['del_end'] - d['del_start']} bp)")
    f, rwt, rko = d["F"], d["R_wt"], d["R_ko"]
    print(f"F     (common)      : {f['seq']:<30} pos {f['start']:>6}-{f['end']:<6} "
          f"len={f['len']:<3} Tm={f['tm']:5.2f} GC={f['gc']:5.1f}% "
          f"self_any_TH={f['self_any']:5.1f} self_end_TH={f['self_end']:5.1f} hairpin_TH={f['hairpin']:5.1f}  "
          f"{_spec_str(f)}")
    print(f"R_wt  (WT-specific) : {rwt['seq']:<30} binds {rwt['start']:>6}-{rwt['end']:<6} "
          f"len={rwt['len']:<3} Tm={rwt['tm']:5.2f} GC={rwt['gc']:5.1f}% "
          f"self_any_TH={rwt['self_any']:5.1f} self_end_TH={rwt['self_end']:5.1f} hairpin_TH={rwt['hairpin']:5.1f}  "
          f"{_spec_str(rwt)}")
    print(f"R_ko  (KO-specific) : {rko['seq']:<30} binds {rko['start']:>6}-{rko['end']:<6} "
          f"len={rko['len']:<3} Tm={rko['tm']:5.2f} GC={rko['gc']:5.1f}% "
          f"self_any_TH={rko['self_any']:5.1f} self_end_TH={rko['self_end']:5.1f} hairpin_TH={rko['hairpin']:5.1f}  "
          f"{_spec_str(rko)}")
    print(f"Product WT (F+R_wt) : {d['product_wt']} bp")
    print(f"Product KO (F+R_ko) : {d['product_ko']} bp")
    print(f"Size difference     : {d['size_diff']} bp")
    print(f"Tm spread (3 primer): {d['tm_spread']:.2f} C")
    print(f"Pair compl ANY_TH   : {d['pair_compl_any']}")
    print(f"Pair compl END_TH   : {d['pair_compl_end']}")
    print()


# ---------------------------------------------------------------------------
# Knock-in (HDR, large cassette insertion) genotyping
# ---------------------------------------------------------------------------
#
# Duality with the deletion/knockout design above: in the KO design, the
# "deletion region" is a span that's present in full_seq (the WT allele) and
# absent from the edited allele. For a knock-in, it's the other way around:
# the cassette is a span that's present in the edited (KI) allele and absent
# from the WT allele. Mathematically this is the exact mirror image, so this
# wrapper builds the KI allele's sequence (WT flanks + cassette spliced in)
# and re-uses design_3primer_genotyping() directly, treating the cassette as
# the "deletion region" — then relabels the output so it makes sense for a
# knock-in:
#     engine's R_wt  (picked *inside* the marked region)      -> our R_ki  (cassette-specific)
#     engine's R_ko  (picked *downstream* of the marked region) -> our R_wt  (WT-specific)
#     engine's product_wt -> our product_ki
#     engine's product_ko -> our product_wt
#
# Three inputs describe the edit: unedited_seq (the WT locus), replaced_seq
# (the part of unedited_seq the cassette replaces — or "NA"/None plus a '^'
# marker in unedited_seq for a pure insertion with nothing removed), and
# cassette_seq (the inserted sequence). resolve_knockin_insertion() turns the
# first two into 0-based coordinates.
#
# Simplest path (recommended): use expand_flanks_from_genome_using_arms()
# with your left/right HDR homology arms instead of figuring out
# unedited_seq/replaced_seq by hand — it locates both arms in the reference
# genome and returns unedited_seq + insert_start/insert_end already resolved
# (whatever sits between the arms in the genome, if anything, is picked up
# automatically), which you then pass straight into
# design_3primer_knockin_genotyping() below via insert_start/insert_end.

def design_3primer_knockin_genotyping(
    unedited_seq,
    replaced_seq,
    cassette_seq,
    insert_start=None,
    insert_end=None,
    product_min=200,
    product_max=3000,
    min_size_diff=500,
    primer_len_range=(18, 28),
    tm_range=(59, 63),
    tm_opt=60.0,
    gc_range=(40, 60),
    gc_opt=50.0,
    max_search_upstream=1600,
    max_search_into_cassette=1600,
    max_search_downstream=1600,
    max_self_any_th=47.0,
    max_self_end_th=47.0,
    max_hairpin_th=47.0,
    max_pair_compl_any_th=47.0,
    max_pair_compl_end_th=47.0,
    max_tm_spread=3.0,
    top_n=5,
    candidate_pool_size=150,
    genome_build=None,
    genome_fastas=None,
    max_genome_hits=1,
    specificity_min_pident=90.0,
    specificity_min_coverage=0.9,
    blastn_bin="blastn",
    makeblastdb_bin="makeblastdb",
):
    """Design a 3-primer genotyping set for an HDR knock-in that inserts a
    (typically large) cassette.

    Args:
        unedited_seq: the *unedited* WT genomic sequence spanning the
            insertion site, ideally with generous flanking sequence already
            included (see expand_flanks_from_genome to auto-fetch this).
        replaced_seq: which part of unedited_seq (if any) gets replaced/removed
            by the cassette.
                - An exact substring of unedited_seq (e.g. a STOP codon being
                  removed while inserting a tag) — its position is found
                  automatically.
                - None / "" / "NA" for a *pure insertion* with nothing
                  removed. In that case mark the exact insertion point inside
                  unedited_seq with a single '^' character (e.g.
                  "...ATGAAG^CTGGAT..."); it's stripped out automatically.
            Ignored if insert_start/insert_end are both given directly (used
            internally after genome flank expansion — see expand_flanks_from_genome).
        cassette_seq: the sequence that gets inserted by HDR (reporter, tag,
            selection cassette, etc.). This does not need to exist in the
            reference genome, and is not genome-BLASTed as a whole — only the
            individual candidate primers picked inside it are.
        insert_start, insert_end: optional pre-resolved 0-based coordinates in
            unedited_seq (skips replaced_seq/marker resolution). Mainly used
            internally when unedited_seq has already been through
            expand_flanks_from_genome.
        (all other args mirror design_3primer_genotyping)

    Returns:
        A list of dicts (up to top_n), sorted ascending by score. Each dict has:
            F               : common forward primer (genomic, upstream of the insertion site)
            R_wt            : WT-specific reverse primer (genomic, downstream of the
                               insertion site; F+R_wt only falls in the amplifiable size
                               range in the WT allele — in the KI allele the cassette
                               pushes this primer too far away)
            R_ki            : KI-specific reverse primer (inside cassette_seq; this
                               binding site doesn't exist at all in the WT allele)
            product_wt, product_ki, size_diff, tm_spread, pair_compl_any, pair_compl_end,
            cassette_start, cassette_end, score
    """
    cassette_seq = cassette_seq.strip()

    if insert_start is not None and insert_end is not None:
        clean_seq = unedited_seq.strip()
    else:
        clean_seq, insert_start, insert_end = resolve_knockin_insertion(unedited_seq, replaced_seq)

    full_seq_ki = clean_seq[:insert_start] + cassette_seq + clean_seq[insert_end:]
    cassette_start = insert_start
    cassette_end = insert_start + len(cassette_seq)

    print(f"[knock-in mode] Internally reusing the deletion/knockout engine: its 'R_wt' "
          f"is your R_ki (cassette-specific), its 'R_ko' is your R_wt (WT-specific), and "
          f"its 'deletion region' below refers to the cassette "
          f"({cassette_start}-{cassette_end}, {len(cassette_seq):,}bp). This only matters "
          f"if you see diagnostic output below.")

    raw_results = design_3primer_genotyping(
        full_seq_ki,
        del_start=cassette_start,
        del_end=cassette_end,
        product_min=product_min,
        product_max=product_max,
        min_size_diff=min_size_diff,
        primer_len_range=primer_len_range,
        tm_range=tm_range,
        tm_opt=tm_opt,
        gc_range=gc_range,
        gc_opt=gc_opt,
        max_search_upstream=max_search_upstream,
        max_search_downstream=max_search_downstream,
        max_search_into_deletion=max_search_into_cassette,
        max_self_any_th=max_self_any_th,
        max_self_end_th=max_self_end_th,
        max_hairpin_th=max_hairpin_th,
        max_pair_compl_any_th=max_pair_compl_any_th,
        max_pair_compl_end_th=max_pair_compl_end_th,
        max_tm_spread=max_tm_spread,
        top_n=top_n,
        candidate_pool_size=candidate_pool_size,
        genome_build=genome_build,
        genome_fastas=genome_fastas,
        max_genome_hits=max_genome_hits,
        specificity_min_pident=specificity_min_pident,
        specificity_min_coverage=specificity_min_coverage,
        blastn_bin=blastn_bin,
        makeblastdb_bin=makeblastdb_bin,
    )

    results = []
    for r in raw_results:
        results.append(dict(
            F=r["F"],
            R_ki=r["R_wt"],
            R_wt=r["R_ko"],
            product_ki=r["product_wt"],
            product_wt=r["product_ko"],
            size_diff=r["size_diff"],
            tm_spread=r["tm_spread"],
            pair_compl_any=dict(F_Rki=r["pair_compl_any"]["F_Rwt"],
                                 F_Rwt=r["pair_compl_any"]["F_Rko"],
                                 Rki_Rwt=r["pair_compl_any"]["Rwt_Rko"]),
            pair_compl_end=dict(F_Rki=r["pair_compl_end"]["F_Rwt"],
                                 F_Rwt=r["pair_compl_end"]["F_Rko"],
                                 Rki_Rwt=r["pair_compl_end"]["Rwt_Rko"]),
            cassette_start=cassette_start,
            cassette_end=cassette_end,
            score=r["score"],
        ))
    return results


def print_knockin_design(d, rank=None):
    """Like print_design(), but labeled for a knock-in (R_ki / product_ki
    instead of R_wt-inside-deletion / product from the deletion design)."""
    print(f"=== Knock-in Design{' #' + str(rank) if rank else ''} (score={d['score']:.2f}) ===")
    print(f"Cassette coordinates (0-based, python slice, within the KI allele sequence): "
          f"{d['cassette_start']}-{d['cassette_end']} "
          f"(length {d['cassette_end'] - d['cassette_start']} bp)")
    f, rwt, rki = d["F"], d["R_wt"], d["R_ki"]
    print(f"F     (common)      : {f['seq']:<30} pos {f['start']:>6}-{f['end']:<6} "
          f"len={f['len']:<3} Tm={f['tm']:5.2f} GC={f['gc']:5.1f}% "
          f"self_any_TH={f['self_any']:5.1f} self_end_TH={f['self_end']:5.1f} hairpin_TH={f['hairpin']:5.1f}  "
          f"{_spec_str(f)}")
    print(f"R_wt  (WT-specific) : {rwt['seq']:<30} binds {rwt['start']:>6}-{rwt['end']:<6} "
          f"len={rwt['len']:<3} Tm={rwt['tm']:5.2f} GC={rwt['gc']:5.1f}% "
          f"self_any_TH={rwt['self_any']:5.1f} self_end_TH={rwt['self_end']:5.1f} hairpin_TH={rwt['hairpin']:5.1f}  "
          f"{_spec_str(rwt)}")
    print(f"R_ki  (KI-specific) : {rki['seq']:<30} binds {rki['start']:>6}-{rki['end']:<6} "
          f"len={rki['len']:<3} Tm={rki['tm']:5.2f} GC={rki['gc']:5.1f}% "
          f"self_any_TH={rki['self_any']:5.1f} self_end_TH={rki['self_end']:5.1f} hairpin_TH={rki['hairpin']:5.1f}  "
          f"{_spec_str(rki)}")
    print(f"Product WT (F+R_wt) : {d['product_wt']} bp")
    print(f"Product KI (F+R_ki) : {d['product_ki']} bp")
    print(f"Size difference     : {d['size_diff']} bp")
    print(f"Tm spread (3 primer): {d['tm_spread']:.2f} C")
    print(f"Pair compl ANY_TH   : {d['pair_compl_any']}")
    print(f"Pair compl END_TH   : {d['pair_compl_end']}")
    print()


# ---------------------------------------------------------------------------
# Usage examples
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Usage 1) supply deletion_seq as a sequence directly:
    #   design_3primer_genotyping(full_seq, deletion_seq=deletion_seq)
    # Usage 2) supply coordinates directly if you already know them:
    #   design_3primer_genotyping(full_seq, del_start=1234, del_end=5678)

    import random
    random.seed(42)

    def _rand_seq(n, gc=0.45):
        out = []
        for _ in range(n):
            out.append(random.choice("GC") if random.random() < gc else random.choice("AT"))
        return "".join(out)

    upstream = _rand_seq(800)
    deletion = _rand_seq(1200)
    downstream = _rand_seq(800)
    full_seq_demo = upstream + deletion + downstream

    print("### Demo: sanity test with a randomly generated sequence (real primer3-py calls) ###\n")
    designs = design_3primer_genotyping(
        full_seq_demo,
        deletion_seq=deletion,
        product_min=200,
        product_max=1500,
        min_size_diff=500,
        top_n=3,
        genome_fastas=[],  # demo only: test primer3 calculations without genome BLAST
    )

    if not designs:
        print("Could not find a design that satisfies the constraints. Try relaxing them.")
    else:
        for i, d in enumerate(designs, 1):
            print_design(d, rank=i)
