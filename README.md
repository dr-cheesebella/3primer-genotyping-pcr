# 3-Primer Genotyping PCR Primer Designer

Automated design of 3-primer genotyping PCR sets (common forward + a primer
specific to each allele) for validating CRISPR/Cas9 edits by gel band size —
covers both a deletion/knockout (KO) and an HDR knock-in (KI, e.g. inserting a
reporter/tag/selection cassette), since the two are the exact mirror image of
each other under the hood.

## How it works

Every edit (a deletion or an insertion) is a "marked region" present on one
allele and absent on the other. Three primers tell the alleles apart by PCR
product size:

- **F** — a common forward primer upstream of the marked region (present on both alleles)
- a primer picked *inside* the marked region — only binds on the allele that has it
- a primer picked *outside/downstream* of the marked region — binds on both alleles' underlying genomic sequence, but only falls in the amplifiable size range on the allele where the marked region is absent (on the other allele, it spans the whole region and is too large to amplify normally)

For a **deletion/knockout (KO)**: the marked region is the deleted sequence
(present in WT, absent in KO). Concretely, F + **R_wt** (inside the deletion,
WT-specific) + **R_ko** (downstream, KO-specific).

For an **HDR knock-in (KI)**: the marked region is the inserted cassette
(present in KI, absent in WT) — the mirror image. Concretely, F + **R_ki**
(inside the cassette, KI-specific) + **R_wt** (downstream, WT-specific).

Either way, it picks a combination where both alleles' PCR products fall
within a target size range and differ enough in size to be distinguished
clearly on a gel.

Optionally, it can also:
- pull extra flanking sequence (default 50 kb each side) around your deletion directly from a reference genome, so the KO-specific primer search isn't limited to whatever sequence you happened to paste in
- check every candidate primer for genome-wide off-target binding sites using local BLAST, and print diagnostics showing exactly why no valid design was found (if that happens)

## What it's built on

- [Primer3](https://github.com/primer3-org/primer3) — via the [`primer3-py`](https://github.com/libnano/primer3-py) Python bindings, for all Tm, GC%, self-complementarity, hairpin, and primer-dimer calculations (the same thermodynamic engine used by Primer3Plus and NCBI Primer-BLAST). No custom approximations are used.
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) (`blastn`, `makeblastdb`, `blastdbcmd`) — for optional genome-wide specificity checks and reference-genome flank extraction, run entirely locally against a reference FASTA you provide (human GRCh37/hg19 or GRCh38/hg38, or mouse GRCm38/mm10 or GRCm39/mm39).

## Setup

Requires Python 3.9 or newer. If your system Python is older, install a newer
one first (e.g. `brew install python` on macOS), then use that to create the
virtual environment below.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install primer3-py
```

For the optional genome-based features (off-target checking, automatic flank
expansion), also install BLAST+ (e.g. `brew install blast` on macOS) and place
a reference genome FASTA under `genomes/` (see `GENOME_BUILD_PATHS` in
`design_3primer_genotyping.py` for expected filenames, or point
`genome_fastas` at your own file).

### About the BLAST index (first run only)

The first time you run the tool against a given genome FASTA, it automatically
builds a BLAST index for it (`makeblastdb`), producing extra files right next
to the FASTA (`.nin`, `.nhr`, `.nsq`, `.ndb`, etc.). This can take a few
minutes for a full human genome, but only happens once — later runs detect and
reuse the existing index.

These index files (and the genome FASTA itself) are intentionally excluded via
`.gitignore` since they're large, auto-generated, and shouldn't be
redistributed. If you clone this repo fresh, you'll need to supply your own
genome FASTA under `genomes/`; the index will rebuild automatically on first
use.

**To set this up from scratch:**

```bash
mkdir -p genomes
cd genomes

# GRCh38 / hg38 (human)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# GRCh37 / hg19 (human, optional, only if you also need this build)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.14_GRCh37.p13/GCA_000001405.14_GRCh37.p13_genomic.fna.gz
gunzip GCA_000001405.14_GRCh37.p13_genomic.fna.gz

# GRCm39 / mm39 (mouse, current/latest assembly -- optional, only if you work with mouse)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
gunzip GCF_000001635.27_GRCm39_genomic.fna.gz

# GRCm38.p6 / mm10 (mouse, older assembly -- optional, kept alongside mm39 since a lot
# of existing mouse data/annotations still use it, similar to hg19 vs hg38)
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
gunzip GCF_000001635.26_GRCm38.p6_genomic.fna.gz

cd ..
```

Only download the builds you actually need — each is a multi-GB file.

Your `genomes/` folder is ready. The BLAST index still needs to be built at
some point (`makeblastdb`) — you just don't have to run that command
**yourself, right now**. Just run the tool as usual (e.g.
`python run_design_KO.py`, `python run_design_KI.py`), and the first time it
touches a given FASTA, it will notice no index exists yet and build one
automatically before continuing (this is the part that can take a few
minutes for a full human genome).

Optionally, if you'd rather build the index ahead of time yourself instead of
waiting for it during your first real run (e.g. kick it off overnight), you
can run the exact same command manually:

```bash
makeblastdb -in genomes/GCF_000001405.40_GRCh38.p14_genomic.fna -dbtype nucl -parse_seqids
makeblastdb -in genomes/GCA_000001405.14_GRCh37.p13_genomic.fna -dbtype nucl -parse_seqids
makeblastdb -in genomes/GCF_000001635.27_GRCm39_genomic.fna -dbtype nucl -parse_seqids
makeblastdb -in genomes/GCF_000001635.26_GRCm38.p6_genomic.fna -dbtype nucl -parse_seqids
```

## Usage

Edit `run_design_KO.py` with your own `FULL_SEQ` and `DELETION_SEQ`, then:

```bash
python run_design_KO.py
```

Or call the function directly:

```python
from design_3primer_genotyping import design_3primer_genotyping, print_design

designs = design_3primer_genotyping(
    full_seq,
    deletion_seq=deletion_seq,
    product_min=200,
    product_max=1500,
    min_size_diff=500,
    genome_build="hg38",   # or "hg19", "both" (human); "mm10"/"mm39" (mouse); pass genome_fastas=[] to skip off-target checking
)
for i, d in enumerate(designs, 1):
    print_design(d, rank=i)
```

### HDR knock-in (large cassette insertion)

For a knock-in that inserts a large cassette (reporter, tag, selection
marker, etc.) instead of a deletion, use `run_design_KI.py`. It takes three
sequences:

- `LEFT_ARM` / `RIGHT_ARM` — your HDR donor's homology arms, in plus-strand
  (reference) orientation. You don't need to know or specify what (if
  anything) sits between them in the genome — the tool locates both arms in
  the reference genome and figures that out automatically (nothing, for a
  pure insertion where the arms are directly adjacent; or a real stretch, if
  your design replaces/removes a short region, e.g. a STOP codon).
- `CASSETTE_SEQ` — the sequence HDR inserts. Doesn't need to exist in the
  reference genome.

```bash
python run_design_KI.py
```

Or call the functions directly:

```python
from design_3primer_genotyping import (
    expand_flanks_from_genome_using_arms,
    design_3primer_knockin_genotyping,
    print_knockin_design,
)

wt_locus_seq, insert_start, insert_end, info = expand_flanks_from_genome_using_arms(
    left_arm, right_arm, genome_build="hg38", flank_kb=50,
)

designs = design_3primer_knockin_genotyping(
    wt_locus_seq,
    replaced_seq=None,   # already resolved via the arms above
    cassette_seq=cassette_seq,
    insert_start=insert_start,
    insert_end=insert_end,
    product_min=200,
    product_max=3000,
    min_size_diff=500,
    genome_build="hg38",
)
for i, d in enumerate(designs, 1):
    print_knockin_design(d, rank=i)
```

(If you'd rather skip the genome lookup and hand-assemble your own WT
sequence, `design_3primer_knockin_genotyping()` also accepts
`unedited_seq` + `replaced_seq` directly — pass an exact substring to
replace, or `"NA"` plus a single `^` character marking the insertion point
in `unedited_seq` for a pure insertion.)

### Running without a reference genome

Both `ref_genome` (in `run_design_KO.py`) and `ref_genome` (in
`run_design_KI.py`) can be left empty (`""` or `None`). This skips *both*
automatic flank expansion *and* genome-wide off-target checking entirely — no
BLAST step runs at all, and the console prints a message confirming this.
Useful if you don't have BLAST+/a genome FASTA set up, or just want a quick
check using the sequence exactly as pasted.

Caveats when running this way:
- Make sure your pasted sequence(s) already include enough flanking sequence — there's no auto-expansion to fall back on.
- For KI, `LEFT_ARM`/`RIGHT_ARM` are assumed to sit directly adjacent (a pure insertion with nothing between them), since there's no genome lookup to determine any gap automatically.

## Author

Sim Sakong, Hansen Lab, MIT — v1, September 2026

## License

MIT (see `LICENSE`).
