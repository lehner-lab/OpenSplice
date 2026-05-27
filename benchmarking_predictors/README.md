# open_splice_benchmarking

Benchmarking scripts for splicing predictors on OpenSPlice.

This repository contains the workflow used for the **Figure 4** section of the OpenSPlice manuscript, including inference, post-processing, and figure recreation.

Associated data: Figshare DOI **10.6084/m9.figshare.32337414**.

## Repository layout
- `scripts/inference/`: inference entrypoints with standardized `<model>_<mode>_<variant>_inference.py` naming.
- `scripts/processing/`: post-processing, aggregation, and merge scripts.
- `scripts/setup_*.sh`: environment bootstrap scripts by model family.
- `scripts/utils/`: one-off helper utilities.
- `docs/metadata/input_column_requirements.tsv`: required input columns by file and consumer script.
- `docs/metadata/inference_input_requirements.md`: human-readable input requirements guidance.
- `docs/requirements/requirements_*.txt`: model-specific Python dependencies.
- `docs/requirements/alphagenome_local_setup.txt`: AlphaGenome local environment notes.
- `plotting/figure_4_and_ED_figures_manuscript_plotting.ipynb`: notebook for recreating manuscript plots.

For additional operational details (naming conventions, setup, run order, output contracts, and release notes), see `docs/metadata/repository_notes.md`.
