#!/usr/bin/env python3
"""
run_design_KO.py  [KO / deletion-knockout genotyping]

Version: v1
Author: Sim Sakong, Hansen Lab, MIT
Date: 2026-09-04

Runner script for calling design_3primer_genotyping() on a real sequence, for
genotyping a CRISPR deletion (knockout) allele.
Fill in the FULL_SEQ / DELETION_SEQ variables below with your actual sequences,
then run:

    python run_design_KO.py

For an HDR knock-in (large cassette insertion) instead, see run_design_KI.py.
"""

from design_3primer_genotyping import (
    design_3primer_genotyping,
    expand_flanks_from_genome,
    print_design,
)

# ---------------------------------------------------------------------------
# Fill in your actual sequences here
# ---------------------------------------------------------------------------

# 1) The full genomic sequence (deletion region + enough flanking sequence on
#    both sides).
#    Note: to design a KO-specific reverse primer you need at least a few
#    hundred bp of sequence past the end of the deletion (so that F+R_ko is
#    too large to amplify in the WT allele).
FULL_SEQ = """
PASTE_YOUR_FULL_SEQUENCE_HERE
"""

# 2) The sequence to be deleted. Must be an exact substring of FULL_SEQ
#    (case and every base must match exactly so find() can locate it).
DELETION_SEQ = """
PASTE_YOUR_DELETION_SEQUENCE_HERE
"""

# 3) Match this to whichever reference you're using: "hg19" / "hg38" / "both"
#    (Note: expand_flanks_from_genome doesn't support "both" — you must pick one.
#     design_3primer_genotyping's genome_build can be "both".)
#    Leave this empty ("" or None) to skip BOTH automatic flank expansion AND
#    genome-wide off-target checking entirely -- the sequence is used exactly
#    as pasted above, with no BLAST step at all.
ref_genome = "hg38"

# 4) How much flanking sequence to automatically pull from the genome on each
#    side of the deletion (in kb) -- only used if ref_genome is set above.
FLANK_KB = 50

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    full_seq = "".join(FULL_SEQ.split())        # strip whitespace/newlines
    deletion_seq = "".join(DELETION_SEQ.split())

    if ref_genome:
        # Step 1: re-fetch sequence from the genome with FLANK_KB of generous
        # flanking on each side of the deletion.
        # (This fixes the KO-specific reverse primer search being too tight —
        #  instead of only the sequence right after the deletion, candidates can
        #  now be chosen anywhere within a 50kb window.)
        expanded_full_seq, del_start, del_end, info = expand_flanks_from_genome(
            full_seq,
            deletion_seq=deletion_seq,
            genome_build=ref_genome,
            flank_kb=FLANK_KB,
        )
        print(f"[run_design_KO] Obtained genome context: {info['chrom']} "
              f"upstream={info['upstream_flank_bp']:,}bp, "
              f"downstream={info['downstream_flank_bp']:,}bp\n")
    else:
        print("[run_design_KO] No reference genome specified (ref_genome is empty) -- "
              "skipping automatic flank expansion and genome-wide off-target checking. "
              "Using the sequence exactly as pasted above; make sure FULL_SEQ already "
              "has enough flanking sequence on both sides.\n")
        expanded_full_seq = full_seq
        del_start = del_end = None   # resolved from deletion_seq inside design_3primer_genotyping instead

    # Step 2: run the 3-primer design on the expanded sequence/coordinates
    designs = design_3primer_genotyping(
        expanded_full_seq,
        deletion_seq=deletion_seq if del_start is None else None,
        del_start=del_start,
        del_end=del_end,
        product_min=200,
        product_max=2000,
        min_size_diff=500,
        genome_build=ref_genome if ref_genome else None,
        genome_fastas=[] if not ref_genome else None,
        top_n=5,
    )

    if not designs:
        print("Could not find a design that satisfies the constraints. "
              "See the diagnostics above and try relaxing the constraints.")
    else:
        for i, d in enumerate(designs, 1):
            print_design(d, rank=i)
