# Repository operational notes

This document contains operational details moved from the top-level `README.md` to keep the landing page concise.

## Naming convention
- Inference scripts: `<model>_<mode>_<variant>_inference.py`
- Processing scripts: `process_<model>_<mode>_outputs.py`

Legacy suffixes like `_no_corrs` were removed to keep names concise and stage-focused.

## Input files
Primary metadata files expected under `data/input/`:
- `opensplice_predictors_benchmarking_variant_metadata.tsv`
- `opensplice_predictors_benchmarking_exon_metadata.tsv`

See `docs/metadata/input_column_requirements.tsv` for required columns per stage.

### Public-repo recommendation for input-file docs
- Keep the input-file requirement docs in the public repo.
- These files are not just internal notes: they define the schema contract needed to reproduce benchmarking outputs.
- If you want a cleaner public landing page, keep the detailed schema in `docs/metadata/` and leave only a short pointer in this README.
- Do **not** include any private or non-redistributable sample data in the repo; publish only schema/column requirements and synthetic examples.

## Setup
Use one environment per model family.

- SpliceAI: `bash scripts/setup_spliceai.sh`
- Pangolin: `bash scripts/setup_pangolin.sh`
- SpliceTransformer: `bash scripts/setup_splicetransformer.sh`
- AlphaGenome: `bash scripts/setup_alphagenome.sh`

Install Python packages from the corresponding file in `docs/requirements/`.

## Recreating manuscript plots
To recreate plots from the manuscript (including Figure 4 and Extended Data figures), use:

- `plotting/figure_4_and_ED_figures_manuscript_plotting.ipynb`

This notebook is the primary plotting entrypoint for figure recreation.

## Output directory contract
- SpliceAI genome inference outputs:
  - `results/spliceai/genome_mode/snvs`
  - `results/spliceai/genome_mode/deletions`
- Pangolin genome inference outputs:
  - `results/pangolin/genome_mode/snvs`
  - `results/pangolin/genome_mode/deletions`
- Minigene outputs follow the same pattern under `results/<model>/minigene_mode/...`.
- SpliceTransformer genomic inference outputs: `results/splice_transformer/genomic/raw_parquet`
- SpliceTransformer minigene inference outputs: `results/splice_transformer/minigene_mode/inference`
- SpliceTransformer processing defaults read these same directories via `ST_DIR`.

## Recommended run order
1. Run scripts in `scripts/inference/` to generate model outputs.
2. Run scripts in `scripts/processing/` to normalize/merge downstream tables.
3. Join curated outputs for plotting/statistics.
4. Recreate figures with the notebook in `plotting/`.

## Public release checklist (suggested)
- Add a LICENSE file (if not already present in your org/release flow).
- Add a CONTRIBUTING guide and issue templates if you want outside contributions.
- Pin expected Python version(s) for each model environment.
- Add one minimal end-to-end command example (single model, tiny synthetic input).
- Ensure all docs reference scripts that currently exist in `scripts/`.
