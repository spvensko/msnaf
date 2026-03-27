# mSNAF

<p align="center">
  <img src="docs/msnaf-example.png" alt="mSNAF example output" width="50%">
</p>

`msnaf` is a minimal, translation-only rewrite of the SNAF execution path you described.

It takes `counts.original.full.txt`, applies the same junction-vs-control filtering logic that SNAF uses, reconstructs splice-junction sequences against the same reference bundle, and emits only:

- `uid`
- `coord`
- `peptide`
- `coding_sequence`
- `peptide_context`

It does not do any of the following:

- HLA parsing
- MHC binding prediction
- immunogenicity scoring
- pickle serialization
- multiprocessing
- T-antigen/B-antigen reporting layers

## How It Differs From SNAF

SNAF is a large framework with T-cell and B-cell antigen workflows. `msnaf` is intentionally narrower.

SNAF-style flow:

1. Read a pruned or reduced junction matrix.
2. Filter junctions against GTEx and optional controls.
3. Translate candidate events.
4. Predict peptide-MHC binding.
5. Score immunogenicity.
6. Generate frequency tables, burden tables, candidate reports, and viewers.

`msnaf` flow:

1. Read `counts.original.full.txt`.
2. Strip AltAnalyze's trailing coordinate suffix from each row id.
3. Filter junctions against GTEx and optional controls using the same SNAF thresholds.
4. Reconstruct splice peptides and their coding sequences.
5. Write one CSV.

## Input Expectations

`msnaf` expects the same extracted reference bundle SNAF uses, including:

- `Alt91_db/Hs_Ensembl_exon_add_col.txt`
- `Alt91_db/mRNA-ExonIDs.txt`
- `Alt91_db/Hs_gene-seq-2000_flank.fa`
- `Alt91_db/Homo_sapiens.GRCh38.91.gtf`
- `Alt91_db/df_start_codon.txt`
- `controls/GTEx_junction_counts.h5ad`

If present, `controls/tcga_matched_control_junction_count.h5ad` is also loaded by default, matching your existing SNAF process.

## Local Run

Create an environment and install:

```bash
cd ~/dev/msnaf
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Run:

```bash
msnaf \
  --counts /path/to/counts.original.full.txt \
  --refs /path/to/SNAF_refs \
  --output /path/to/snaf_peptides.csv
```

Optional flags:

```bash
--strict
--filter-mode maxmin
--filter-mode prevalance
--not-in-db
--skip-tcga-control
```

## Output

Example CSV rows:

```text
uid,coord,peptide,coding_sequence,peptide_context
ENSG00000110427:E7.1-E8.1,chr11:33559911-33561676(+),PLEYPNLDI,CCACTGGAATATCCCAACCTTGACATAT,PLEYPNLDISETTRDYWVI
ENSG00000110427:E18.1-E20.1,chr11:33583501-33591237(+),KPVQGFDYA,AAGCCTGTGCAAGGCTTTGATTATGCCA,LPQRAKPVQGFDYAKQHLGQQGAD
```

## Docker

Build from the `msnaf` repo root:

```bash
cd ~/dev/msnaf
docker build -t msnaf:latest .
```

Run:

```bash
docker run --rm \
  -v /path/to/data:/data \
  msnaf:latest \
  --counts /data/counts.original.full.txt \
  --refs /data/refs \
  --output /data/snaf_peptides.csv
```

The container entrypoint is `msnaf`, so the runtime arguments are the CLI flags directly.

## Testing

```bash
cd ~/dev/msnaf
python3 -m unittest discover -s tests -v
```
