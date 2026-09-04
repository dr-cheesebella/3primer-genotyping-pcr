#!/usr/bin/env python3
"""
run_design_KI.py  [KI / HDR knock-in genotyping, large cassette insertion]

Version: v1
Author: Sim Sakong, Hansen Lab, MIT
Date: 2026-09-04

Runner script for genotyping an HDR knock-in that inserts a large cassette
(reporter, tag, selection marker, etc.) at a specific site.
Fill in LEFT_ARM / RIGHT_ARM / CASSETTE_SEQ below with your actual sequences,
then run:

    python run_design_KI.py

For a deletion/knockout instead, see run_design.py (a.k.a. run_design_KO.py).
"""

from design_3primer_genotyping import (
    expand_flanks_from_genome_using_arms,
    design_3primer_knockin_genotyping,
    print_knockin_design,
)

# ---------------------------------------------------------------------------
# Fill in your actual sequences here
# ---------------------------------------------------------------------------

# 1) The left and right HDR homology arms, exactly as used in your donor
#    construct (plus-strand / reference orientation -- if your design is on
#    the minus strand, reverse-complement both arms first).
#
#    You don't need to know or specify what (if anything) sits between them
#    in the genome -- expand_flanks_from_genome_using_arms() below locates
#    both arms and figures that out automatically (zero-length for a pure
#    insertion where the arms are directly adjacent, or a real stretch if
#    your design replaces/removes a short region, e.g. a STOP codon).
LEFT_ARM = """
PASTE_YOUR_LEFT_HOMOLOGY_ARM_SEQUENCE_HERE
"""

RIGHT_ARM = """
PASTE_YOUR_RIGHT_HOMOLOGY_ARM_SEQUENCE_HERE
"""

# 2) The cassette sequence that HDR inserts (reporter, tag, selection marker,
#    etc.). This does NOT need to exist in the reference genome -- it's your
#    construct's sequence, not a genomic sequence.
CASSETTE_SEQ = """
PASTE_YOUR_CASSETTE_SEQUENCE_HERE
"""

# 3) Match this to whichever reference your arms came from: "hg19" / "hg38"
#    Leave this empty ("" or None) to skip BOTH automatic flank expansion AND
#    genome-wide off-target checking entirely. Without a genome, LEFT_ARM and
#    RIGHT_ARM are assumed to sit directly adjacent (a pure insertion with
#    nothing between them) -- if your design actually replaces/removes a
#    stretch of sequence, you'll need to include it yourself (e.g. append it
#    to LEFT_ARM), and both arms should already include enough flanking
#    sequence for primer search room.
ref_genome = "hg38"

# 4) How much flanking sequence to automatically pull from the genome beyond
#    the outer edges of the arms (in kb) -- gives the WT-specific reverse
#    primer more room to be found away from the cassette junction. Only used
#    if ref_genome is set above.
FLANK_KB = 50

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    left_arm = "".join(LEFT_ARM.split())      # strip whitespace/newlines
    right_arm = "".join(RIGHT_ARM.split())
    cassette_seq = "".join(CASSETTE_SEQ.split())

    if ref_genome:
        # Step 1: locate both homology arms in the reference genome and fetch a
        # single contiguous WT-locus sequence with FLANK_KB of generous flanking
        # beyond each arm. Whatever (if anything) sits between the arms in the
        # genome is picked up automatically -- no need to specify it by hand.
        wt_locus_seq, insert_start, insert_end, info = expand_flanks_from_genome_using_arms(
            left_arm,
            right_arm,
            genome_build=ref_genome,
            flank_kb=FLANK_KB,
        )
        print(f"[run_design_KI] Obtained genome context: {info['chrom']} "
              f"left_arm={info['left_arm_genomic']}, right_arm={info['right_arm_genomic']}, "
              f"gap_between_arms={info['gap_len']:,}bp, "
              f"upstream_flank={info['upstream_flank_bp']:,}bp, "
              f"downstream_flank={info['downstream_flank_bp']:,}bp\n")
    else:
        print("[run_design_KI] No reference genome specified (ref_genome is empty) -- "
              "skipping automatic flank expansion and genome-wide off-target checking. "
              "LEFT_ARM and RIGHT_ARM are assumed to sit directly adjacent (a pure "
              "insertion with nothing between them); make sure both arms already "
              "include enough flanking sequence.\n")
        wt_locus_seq = left_arm + right_arm
        insert_start = insert_end = len(left_arm)

    # Step 2: run the 3-primer knock-in design on the expanded WT sequence,
    # splicing in the cassette at the (now resolved) insertion coordinates.
    designs = design_3primer_knockin_genotyping(
        wt_locus_seq,
        replaced_seq=None,          # already resolved above
        cassette_seq=cassette_seq,
        insert_start=insert_start,
        insert_end=insert_end,
        product_min=200,
        product_max=3000,
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
            print_knockin_design(d, rank=i)
