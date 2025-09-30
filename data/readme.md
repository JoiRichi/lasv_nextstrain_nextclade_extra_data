Data dictionary and provenance for the `data/` directory.

## mammarena_test/

- Contents: Mammarenavirus sequence sets (NCBI TaxID 1653394) prepared for testing.
- Composition: GenBank accessions only (no RefSeq).
- Date range: 25 Sep 1900 – 27 Sep 2025.
- Length filter: minimum sequence length 600 nt.
- Segments: S and L provided separately as FASTA files.

Files
- `all_mammarena_S_test.fasta` — S segment sequences.
- `all_mammarena_L_test.fasta` — L segment sequences.

## nextclade_output/

CSV outputs from Nextclade runs using the above test inputs. Each file follows standard Nextclade CSV schema (metrics, QC flags, gene‑wise results where applicable).

Files
- `nextclade_S_segment_mammarena_test.csv`
- `nextclade_L_segment_mammarena_test.csv`
- `nextclade_GPC_mammarena_test.csv`

## nextstrain_trees_30_09_25/

JSON exports of LASV builds compatible with Auspice visualization.

Files
- `lassa_s.json` — LASV S segment tree.
- `lassa_l.json` — LASV L segment tree.
- `lassa_gpc.json` — LASV GPC tree.
- `nextstrain_meta_csv_30_09_25/` — associated metadata CSVs for the above exports.

## ncbi_virus/

Snapshot of LASV metadata downloaded from NCBI Virus (as of 30/09/2025).

Files
- `all_lasv_ncbi_virus_metadata_30_09_25.csv`

## lassa-main/

Read‑only copy of `https://github.com/nextstrain/lassa`, included to ensure long‑term reproducibility of manuscript results. Do not modify.