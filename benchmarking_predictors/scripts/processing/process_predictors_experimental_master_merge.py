"""
Final merge of archived-ID predictor outputs into the OpenSplice merged table.

Recommended input:
- Supplementary Table 12 from the manuscript repository:
  `Supplementary_Table_12_OpenSplice_dataset_with_splicing_variant_effect_predictions.tsv`

This table contains:
- variant metadata
- experimental PSI values
- model predictions at canonical splice sites

for 590,104 variants across 599 exons where all models successfully produced
predictions.

Before running:
1. Download Supplementary Table 12 from the manuscript repository.
2. Provide your local path via `--experimental-master-file`.
3. Choose your output path with `--out-file`.

Important:
If the input table already contains predictor score/output columns from a previous
merge, this script automatically removes those columns first, then re-merges all
predictor files to regenerate them cleanly.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


# ============================================================
# DEFAULT INPUT / OUTPUT PATHS
# ============================================================

DEFAULT_EXPERIMENTAL_MASTER_FILE = Path("Supplementary_Table_12_OpenSplice_dataset_with_splicing_variant_effect_predictions.tsv")
DEFAULT_OUT_FILE_SUFFIX = "_with_predictors.tsv"


# ============================================================
# DEFAULT PREDICTOR INPUT FILES
# ============================================================

DEFAULT_PREDICTOR_FILES = {
    "SAI_MINIGENE": Path("df_spliceai_minigene_clean_curated.csv"),
    "SAI_GENOME": Path("df_spliceai_genome_canonical_clean_curated.csv"),
    "PANGO_MINIGENE": Path("df_pangolin_minigene_clean_curated.csv"),
    "PANGO_GENOME": Path("df_pangolin_genome_clean_curated.csv"),
    "ST_MINI_SNV": Path("splicetransformer_minigene_snvs_clean_curated.tsv"),
    "ST_MINI_DEL": Path("splicetransformer_minigene_deletions_clean_curated.tsv"),
    "ST_GEN_SNV": Path("splicetransformer_genome_snvs_clean_curated.tsv"),
    "ST_GEN_DEL": Path("splicetransformer_genome_deletions_clean_curated.tsv"),
    "AG_GEN_SNV": Path("alphagenome_genome_snvs_clean_curated.csv"),
    "AG_GEN_DEL": Path("alphagenome_genome_deletions_clean_curated.csv"),
    "AG_MINIGENE": Path("alphagenome_minigene_clean_curated.tsv"),
}


# ============================================================
# COLUMNS TO KEEP FROM EACH PREDICTOR FILE
# ============================================================

SELECTED_COLS = {
    "SAI_MINIGENE": [
        "ensembl_exon_id",
        "variant_id",
        "spliceai_minigene_var_source",
        "spliceai_minigene_acceptor_wt",
        "spliceai_minigene_donor_wt",
        "spliceai_minigene_acceptor_mut",
        "spliceai_minigene_donor_mut",
        "spliceai_minigene_delta_acceptor",
        "spliceai_minigene_delta_donor",
        "spliceai_minigene_delta_mean",
    ],
    "SAI_GENOME": [
        "ensembl_exon_id",
        "variant_id",
        "spliceai_genome_var_source",
        "spliceai_genome_acceptor_wt",
        "spliceai_genome_donor_wt",
        "spliceai_genome_acceptor_mut",
        "spliceai_genome_donor_mut",
        "spliceai_genome_delta_acceptor",
        "spliceai_genome_delta_donor",
        "spliceai_genome_delta_mean",
    ],
    "PANGO_MINIGENE": [
        "ensembl_exon_id",
        "variant_id",
        "pangolin_minigene_var_source",
        "pangolin_minigene_max_signed_delta_acceptor",
        "pangolin_minigene_max_signed_delta_donor",
        "pangolin_minigene_mean_delta_signed",
        "pangolin_minigene_brain_acceptor_wt",
        "pangolin_minigene_brain_donor_wt",
        "pangolin_minigene_brain_acceptor_mut",
        "pangolin_minigene_brain_donor_mut",
        "pangolin_minigene_brain_delta_acceptor",
        "pangolin_minigene_brain_delta_donor",
        "pangolin_minigene_brain_delta_mean",
        "pangolin_minigene_heart_acceptor_wt",
        "pangolin_minigene_heart_donor_wt",
        "pangolin_minigene_heart_acceptor_mut",
        "pangolin_minigene_heart_donor_mut",
        "pangolin_minigene_heart_delta_acceptor",
        "pangolin_minigene_heart_delta_donor",
        "pangolin_minigene_heart_delta_mean",
        "pangolin_minigene_liver_acceptor_wt",
        "pangolin_minigene_liver_donor_wt",
        "pangolin_minigene_liver_acceptor_mut",
        "pangolin_minigene_liver_donor_mut",
        "pangolin_minigene_liver_delta_acceptor",
        "pangolin_minigene_liver_delta_donor",
        "pangolin_minigene_liver_delta_mean",
        "pangolin_minigene_testis_acceptor_wt",
        "pangolin_minigene_testis_donor_wt",
        "pangolin_minigene_testis_acceptor_mut",
        "pangolin_minigene_testis_donor_mut",
        "pangolin_minigene_testis_delta_acceptor",
        "pangolin_minigene_testis_delta_donor",
        "pangolin_minigene_testis_delta_mean",
    ],
    "PANGO_GENOME": [
        "ensembl_exon_id",
        "variant_id",
        "pangolin_genome_max_signed_delta_acceptor",
        "pangolin_genome_max_signed_delta_donor",
        "pangolin_genome_mean_delta_signed",
        "pangolin_genome_brain_acceptor_wt",
        "pangolin_genome_brain_donor_wt",
        "pangolin_genome_brain_acceptor_mut",
        "pangolin_genome_brain_donor_mut",
        "pangolin_genome_brain_delta_acceptor",
        "pangolin_genome_brain_delta_donor",
        "pangolin_genome_brain_delta_mean",
        "pangolin_genome_heart_acceptor_wt",
        "pangolin_genome_heart_donor_wt",
        "pangolin_genome_heart_acceptor_mut",
        "pangolin_genome_heart_donor_mut",
        "pangolin_genome_heart_delta_acceptor",
        "pangolin_genome_heart_delta_donor",
        "pangolin_genome_heart_delta_mean",
        "pangolin_genome_liver_acceptor_wt",
        "pangolin_genome_liver_donor_wt",
        "pangolin_genome_liver_acceptor_mut",
        "pangolin_genome_liver_donor_mut",
        "pangolin_genome_liver_delta_acceptor",
        "pangolin_genome_liver_delta_donor",
        "pangolin_genome_liver_delta_mean",
        "pangolin_genome_testis_acceptor_wt",
        "pangolin_genome_testis_donor_wt",
        "pangolin_genome_testis_acceptor_mut",
        "pangolin_genome_testis_donor_mut",
        "pangolin_genome_testis_delta_acceptor",
        "pangolin_genome_testis_delta_donor",
        "pangolin_genome_testis_delta_mean",
    ],
    "ST_MINI": [
        "exon_id",
        "variant_id",
        "splice_transformer_minigene_acceptor_wt",
        "splice_transformer_minigene_donor_wt",
        "splice_transformer_minigene_acceptor_mut",
        "splice_transformer_minigene_donor_mut",
        "splice_transformer_minigene_delta_acceptor",
        "splice_transformer_minigene_delta_donor",
        "splice_transformer_minigene_delta_mean",
        "splice_transformer_minigene_var_source",
    ],
    "ST_GENOME": [
        "exon_id",
        "variant_id",
        "splice_transformer_genomic_acceptor_wt",
        "splice_transformer_genomic_donor_wt",
        "splice_transformer_genomic_acceptor_mut",
        "splice_transformer_genomic_donor_mut",
        "splice_transformer_genomic_delta_acceptor",
        "splice_transformer_genomic_delta_donor",
        "splice_transformer_genomic_delta_mean",
        "splice_transformer_genomic_var_source",
    ],
    "AG_GENOME": [
        "exon_id",
        "variant_id",
        "acceptor_ref_at_canonical",
        "acceptor_alt_at_canonical",
        "donor_ref_at_canonical",
        "donor_alt_at_canonical",
        "delta_acceptor",
        "delta_donor",
        "mean_delta_splice",
    ],
    "AG_MINIGENE": [
        "ensembl_exon_id",
        "Identifier",
        "alphagenome_minigene_acceptor_wt",
        "alphagenome_minigene_donor_wt",
        "alphagenome_minigene_acceptor_mut",
        "alphagenome_minigene_donor_mut",
        "alphagenome_minigene_delta_acceptor",
        "alphagenome_minigene_delta_donor",
        "alphagenome_minigene_delta_mean",
    ],
}

PREFIX = {
    "SAI_MINIGENE": "spliceai_minigene__",
    "SAI_GENOME": "spliceai_genome__",
    "PANGO_MINIGENE": "pangolin_minigene__",
    "PANGO_GENOME": "pangolin_genome__",
    "ST_MINI": "splicetransformer_minigene__",
    "ST_GENOME": "splicetransformer_genome__",
    "AG_GENOME": "alphagenome_genome__",
    "AG_MINIGENE": "alphagenome_minigene__",
}

SCORE_COLS = [
    "spliceai_minigene__spliceai_minigene_delta_mean",
    "spliceai_genome__spliceai_genome_delta_mean",
    "pangolin_minigene__pangolin_minigene_mean_delta_signed",
    "pangolin_genome__pangolin_genome_mean_delta_signed",
    "splicetransformer_minigene__splice_transformer_minigene_delta_mean",
    "splicetransformer_genome__splice_transformer_genomic_delta_mean",
    "alphagenome_genome__mean_delta_splice",
    "alphagenome_minigene__alphagenome_minigene_delta_mean",
]


# ============================================================
# HELPERS
# ============================================================

def load_table(path: Path) -> pd.DataFrame:
    """Load a CSV/TSV table using file extension to infer separator."""
    path = Path(path)

    if path.suffix == ".tsv":
        return pd.read_csv(path, sep="\t", low_memory=False)

    return pd.read_csv(path, low_memory=False)


def extract_variant_token(x: object) -> object:
    """
    Extract the final underscore-delimited token from a variant ID.

    This keeps the merge robust when old/new files use different full variant_id
    prefixes but share the final mutation token.
    """
    if pd.isna(x):
        return pd.NA

    s = str(x).strip()
    return s.rsplit("_", 1)[-1] if "_" in s else s


def prefix_payload_columns(
    df: pd.DataFrame,
    prefix: str,
    protected_cols: Iterable[str],
) -> pd.DataFrame:
    """Prefix all non-key payload columns to avoid collisions after merging."""
    protected = set(protected_cols)

    rename_map = {
        c: f"{prefix}{c}"
        for c in df.columns
        if c not in protected
    }

    return df.rename(columns=rename_map)


def print_duplicate_examples(
    df: pd.DataFrame,
    key_cols: list[str],
    label: str,
    max_rows: int = 10,
) -> None:
    """Print a small duplicate-key preview without relying on notebook display()."""
    dup_preview = df.loc[df.duplicated(key_cols, keep=False), key_cols].head(max_rows)

    print(f"\n{label}")
    print(dup_preview.to_string(index=False))


def prepare_predictor_table(
    path: Path,
    table_name: str,
    selected_cols: list[str],
    prefix: str,
    predictor_id_col: str,
    variant_id_col: str,
) -> pd.DataFrame:
    """
    Prepare one predictor table for archive-aware merging.

    Important:
    - predictor_id_col is the ID column as stored in the old predictor file.
    - It is renamed to `{table_name}__predictor_id` so it cannot collide
      with current master columns such as exon_id or ensembl_exon_id.
    """
    df = load_table(path)

    missing = [c for c in selected_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{table_name} is missing selected columns:\n{missing}\n\n"
            f"Available columns:\n{list(df.columns)}"
        )

    df = df[selected_cols].copy()

    predictor_merge_id = f"{table_name}__predictor_id"

    df[predictor_id_col] = df[predictor_id_col].astype("string")
    df[predictor_merge_id] = df[predictor_id_col]
    df["variant_token"] = df[variant_id_col].map(extract_variant_token).astype("string")

    # Drop original unrenamed predictor ID so it never overwrites master exon_id.
    if predictor_id_col in df.columns:
        df = df.drop(columns=[predictor_id_col])

    key_cols = [predictor_merge_id, "variant_token"]

    n_loaded = len(df)

    df = df.dropna(subset=key_cols).copy()

    n_dup = int(df.duplicated(key_cols).sum())
    if n_dup > 0:
        print(f"\n{table_name}: duplicate predictor keys before dedup: {n_dup}")
        print_duplicate_examples(df, key_cols, label="Duplicate-key preview:")
        df = df.drop_duplicates(key_cols, keep="first").copy()

    protected_cols = [predictor_merge_id, "variant_token"]
    df = prefix_payload_columns(df, prefix=prefix, protected_cols=protected_cols)

    print(f"\n{table_name}")
    print("  loaded rows      :", n_loaded)
    print("  merge-ready rows :", len(df))
    print("  unique keys      :", df[key_cols].drop_duplicates().shape[0])
    print("  predictor ID col :", predictor_id_col)
    print("  renamed key col  :", predictor_merge_id)

    return df


def prepare_concat_predictor_table(
    path_a: Path,
    path_b: Path,
    table_name: str,
    selected_cols: list[str],
    prefix: str,
    predictor_id_col: str,
    variant_id_col: str,
) -> pd.DataFrame:
    """Prepare SNV + deletion predictor files together."""
    df_a = prepare_predictor_table(
        path=path_a,
        table_name=f"{table_name}_SNV",
        selected_cols=selected_cols,
        prefix="",
        predictor_id_col=predictor_id_col,
        variant_id_col=variant_id_col,
    )

    df_b = prepare_predictor_table(
        path=path_b,
        table_name=f"{table_name}_DEL",
        selected_cols=selected_cols,
        prefix="",
        predictor_id_col=predictor_id_col,
        variant_id_col=variant_id_col,
    )

    # Both file-specific keys are normalised to one combined predictor key.
    snv_key = f"{table_name}_SNV__predictor_id"
    del_key = f"{table_name}_DEL__predictor_id"
    combined_key = f"{table_name}__predictor_id"

    df_a = df_a.rename(columns={snv_key: combined_key})
    df_b = df_b.rename(columns={del_key: combined_key})

    key_cols = [combined_key, "variant_token"]

    df = pd.concat([df_a, df_b], ignore_index=True, sort=False)
    n_before = len(df)

    n_dup = int(df.duplicated(key_cols).sum())
    if n_dup > 0:
        print(f"\n{table_name}: duplicate keys after SNV+DEL concat: {n_dup}")
        print_duplicate_examples(df, key_cols, label="Duplicate-key preview:")
        df = df.drop_duplicates(key_cols, keep="first").copy()

    protected_cols = [combined_key, "variant_token"]
    df = prefix_payload_columns(df, prefix=prefix, protected_cols=protected_cols)

    print(f"\n{table_name} combined")
    print("  rows before dedup:", n_before)
    print("  rows after dedup :", len(df))
    print("  unique keys      :", df[key_cols].drop_duplicates().shape[0])
    print("  combined key col :", combined_key)

    return df


def merge_predictor_into_master(
    master: pd.DataFrame,
    pred: pd.DataFrame,
    label: str,
    master_archive_key: str,
    predictor_merge_key: str,
    score_col: str,
) -> pd.DataFrame:
    """
    Left-merge predictor predictions into corrected master.

    Master keeps:
      - exon_id
      - ensembl_exon_id
      - exon_id_archive
      - ensembl_exon_id_archive

    Predictor key is temporary and dropped after merge.
    """
    before_rows = len(master)

    out = master.merge(
        pred,
        left_on=[master_archive_key, "variant_token"],
        right_on=[predictor_merge_key, "variant_token"],
        how="left",
        validate="many_to_one",
    )

    if len(out) != before_rows:
        raise ValueError(f"{label}: row count changed after merge.")

    out = out.drop(columns=[predictor_merge_key])

    n_nonnull = out[score_col].notna().sum() if score_col in out.columns else "score column not found"

    print(f"\nMerged {label}")
    print("  master archive key :", master_archive_key)
    print("  predictor merge key:", predictor_merge_key)
    print("  rows retained      :", len(out))
    print("  non-null score     :", n_nonnull)

    for required_col in ["exon_id", "ensembl_exon_id", "exon_id_archive", "ensembl_exon_id_archive"]:
        if required_col not in out.columns:
            raise KeyError(f"{label}: required master column was lost: {required_col}")

    return out


def load_experimental_master(experimental_master_file: Path) -> pd.DataFrame:
    """Load and validate the corrected experimental master table."""
    final_df = pd.read_csv(experimental_master_file, sep="\t", low_memory=False)

    required_master_cols = [
        "exon_id",
        "ensembl_exon_id",
        "exon_id_archive",
        "ensembl_exon_id_archive",
        "variant_id",
    ]

    missing_master_cols = [c for c in required_master_cols if c not in final_df.columns]
    if missing_master_cols:
        raise ValueError(f"experimental_data_master missing columns: {missing_master_cols}")

    for c in required_master_cols:
        final_df[c] = final_df[c].astype("string")

    if "variant_token" not in final_df.columns:
        final_df["variant_token"] = final_df["variant_id"].map(extract_variant_token).astype("string")
    else:
        final_df["variant_token"] = final_df["variant_token"].astype("string")

    print("=" * 100)
    print("LOADED CORRECTED EXPERIMENTAL MASTER")
    print("=" * 100)
    print("final_df shape:", final_df.shape)
    print("Unique current exon_id:", final_df["exon_id"].nunique())
    print("Unique current ensembl_exon_id:", final_df["ensembl_exon_id"].nunique())
    print("Unique archive exon_id:", final_df["exon_id_archive"].nunique())
    print("Unique archive ensembl_exon_id:", final_df["ensembl_exon_id_archive"].nunique())

    return final_df




def drop_existing_predictor_columns(final_df: pd.DataFrame) -> pd.DataFrame:
    """Drop existing predictor-output columns so re-merging starts from a clean table."""
    predictor_prefixes = tuple(PREFIX.values())
    protected_cols = {
        "exon_id",
        "ensembl_exon_id",
        "exon_id_archive",
        "ensembl_exon_id_archive",
        "variant_id",
        "variant_token",
    }

    cols_to_drop = [
        c
        for c in final_df.columns
        if c.startswith(predictor_prefixes) and c not in protected_cols
    ]

    if cols_to_drop:
        print("\nDetected existing predictor columns in master input; dropping before merge:")
        print(f"  columns to drop: {len(cols_to_drop)}")
        final_df = final_df.drop(columns=cols_to_drop).copy()

    return final_df

def prepare_all_predictors(predictor_files: dict[str, Path]) -> dict[str, pd.DataFrame]:
    """Prepare all predictor tables used in the final merge."""
    return {
        "SAI_MINIGENE": prepare_predictor_table(
            path=predictor_files["SAI_MINIGENE"],
            table_name="SAI_MINIGENE",
            selected_cols=SELECTED_COLS["SAI_MINIGENE"],
            prefix=PREFIX["SAI_MINIGENE"],
            predictor_id_col="ensembl_exon_id",
            variant_id_col="variant_id",
        ),
        "SAI_GENOME": prepare_predictor_table(
            path=predictor_files["SAI_GENOME"],
            table_name="SAI_GENOME",
            selected_cols=SELECTED_COLS["SAI_GENOME"],
            prefix=PREFIX["SAI_GENOME"],
            predictor_id_col="ensembl_exon_id",
            variant_id_col="variant_id",
        ),
        "PANGO_MINIGENE": prepare_predictor_table(
            path=predictor_files["PANGO_MINIGENE"],
            table_name="PANGO_MINIGENE",
            selected_cols=SELECTED_COLS["PANGO_MINIGENE"],
            prefix=PREFIX["PANGO_MINIGENE"],
            predictor_id_col="ensembl_exon_id",
            variant_id_col="variant_id",
        ),
        "PANGO_GENOME": prepare_predictor_table(
            path=predictor_files["PANGO_GENOME"],
            table_name="PANGO_GENOME",
            selected_cols=SELECTED_COLS["PANGO_GENOME"],
            prefix=PREFIX["PANGO_GENOME"],
            predictor_id_col="ensembl_exon_id",
            variant_id_col="variant_id",
        ),
        "ST_MINIGENE": prepare_concat_predictor_table(
            path_a=predictor_files["ST_MINI_SNV"],
            path_b=predictor_files["ST_MINI_DEL"],
            table_name="ST_MINIGENE",
            selected_cols=SELECTED_COLS["ST_MINI"],
            prefix=PREFIX["ST_MINI"],
            predictor_id_col="exon_id",
            variant_id_col="variant_id",
        ),
        "ST_GENOME": prepare_concat_predictor_table(
            path_a=predictor_files["ST_GEN_SNV"],
            path_b=predictor_files["ST_GEN_DEL"],
            table_name="ST_GENOME",
            selected_cols=SELECTED_COLS["ST_GENOME"],
            prefix=PREFIX["ST_GENOME"],
            predictor_id_col="exon_id",
            variant_id_col="variant_id",
        ),
        "AG_GENOME": prepare_concat_predictor_table(
            path_a=predictor_files["AG_GEN_SNV"],
            path_b=predictor_files["AG_GEN_DEL"],
            table_name="AG_GENOME",
            selected_cols=SELECTED_COLS["AG_GENOME"],
            prefix=PREFIX["AG_GENOME"],
            predictor_id_col="exon_id",  # AlphaGenome genome exon_id is Ensembl-like ID.
            variant_id_col="variant_id",
        ),
        "AG_MINIGENE": prepare_predictor_table(
            path=predictor_files["AG_MINIGENE"],
            table_name="AG_MINIGENE",
            selected_cols=SELECTED_COLS["AG_MINIGENE"],
            prefix=PREFIX["AG_MINIGENE"],
            predictor_id_col="ensembl_exon_id",
            variant_id_col="Identifier",
        ),
    }


def merge_all_predictors(
    final_df: pd.DataFrame,
    predictors: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """Merge all prepared predictor tables onto the corrected master table."""
    merge_plan = [
        {
            "pred_key": "SAI_MINIGENE",
            "label": "SAI_MINIGENE",
            "master_archive_key": "ensembl_exon_id_archive",
            "predictor_merge_key": "SAI_MINIGENE__predictor_id",
            "score_col": "spliceai_minigene__spliceai_minigene_delta_mean",
        },
        {
            "pred_key": "SAI_GENOME",
            "label": "SAI_GENOME",
            "master_archive_key": "ensembl_exon_id_archive",
            "predictor_merge_key": "SAI_GENOME__predictor_id",
            "score_col": "spliceai_genome__spliceai_genome_delta_mean",
        },
        {
            "pred_key": "PANGO_MINIGENE",
            "label": "PANGO_MINIGENE",
            "master_archive_key": "ensembl_exon_id_archive",
            "predictor_merge_key": "PANGO_MINIGENE__predictor_id",
            "score_col": "pangolin_minigene__pangolin_minigene_mean_delta_signed",
        },
        {
            "pred_key": "PANGO_GENOME",
            "label": "PANGO_GENOME",
            "master_archive_key": "ensembl_exon_id_archive",
            "predictor_merge_key": "PANGO_GENOME__predictor_id",
            "score_col": "pangolin_genome__pangolin_genome_mean_delta_signed",
        },
        {
            "pred_key": "ST_MINIGENE",
            "label": "ST_MINIGENE",
            "master_archive_key": "exon_id_archive",
            "predictor_merge_key": "ST_MINIGENE__predictor_id",
            "score_col": "splicetransformer_minigene__splice_transformer_minigene_delta_mean",
        },
        {
            "pred_key": "ST_GENOME",
            "label": "ST_GENOME",
            "master_archive_key": "exon_id_archive",
            "predictor_merge_key": "ST_GENOME__predictor_id",
            "score_col": "splicetransformer_genome__splice_transformer_genomic_delta_mean",
        },
        {
            "pred_key": "AG_GENOME",
            "label": "AG_GENOME",
            "master_archive_key": "ensembl_exon_id_archive",
            "predictor_merge_key": "AG_GENOME__predictor_id",
            "score_col": "alphagenome_genome__mean_delta_splice",
        },
        {
            "pred_key": "AG_MINIGENE",
            "label": "AG_MINIGENE",
            "master_archive_key": "ensembl_exon_id_archive",
            "predictor_merge_key": "AG_MINIGENE__predictor_id",
            "score_col": "alphagenome_minigene__alphagenome_minigene_delta_mean",
        },
    ]

    for step in merge_plan:
        final_df = merge_predictor_into_master(
            master=final_df,
            pred=predictors[step["pred_key"]],
            label=step["label"],
            master_archive_key=step["master_archive_key"],
            predictor_merge_key=step["predictor_merge_key"],
            score_col=step["score_col"],
        )

    return final_df


def print_final_summary(final_df: pd.DataFrame) -> None:
    """Print final merged-table dimensions, ID counts, and predictor coverage."""
    print("\n" + "=" * 100)
    print("FINAL MERGED TABLE SUMMARY")
    print("=" * 100)
    print("final_df shape:", final_df.shape)
    print("Unique current exon_id:", final_df["exon_id"].nunique())
    print("Unique current ensembl_exon_id:", final_df["ensembl_exon_id"].nunique())
    print("Unique archive exon_id:", final_df["exon_id_archive"].nunique())
    print("Unique archive ensembl_exon_id:", final_df["ensembl_exon_id_archive"].nunique())

    print("\nCoverage:")
    for c in SCORE_COLS:
        if c in final_df.columns:
            print(f"  {c}: {final_df[c].notna().sum()}")

    available_score_cols = [c for c in SCORE_COLS if c in final_df.columns]

    if available_score_cols:
        rows_with_predictor = int(final_df[available_score_cols].notna().any(axis=1).sum())
    else:
        rows_with_predictor = 0

    print("\nRows with at least one predictor:", rows_with_predictor)

    print("\nFinal columns containing exon_id:")
    for c in final_df.columns:
        if "exon_id" in c:
            print(" ", c)

    print("\nPreview:")
    print(final_df.head().to_string(index=False))


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Final merge of SpliceAI, Pangolin, SpliceTransformer, and AlphaGenome "
            "predictor outputs into the corrected experimental master table."
        )
    )

    parser.add_argument(
        "--experimental-master-file",
        type=Path,
        default=DEFAULT_EXPERIMENTAL_MASTER_FILE,
        help="Path to Supplementary Table 12 (or compatible merged master table) TSV.",
    )
    parser.add_argument(
        "--out-file",
        type=Path,
        default=None,
        help="Path for merged output TSV.",
    )

    parser.add_argument(
        "--predictors-dir",
        type=Path,
        default=None,
        help="Optional base directory for predictor files. Default filenames are resolved under this directory.",
    )

    # Optional overrides for predictor paths.
    for key, default_path in DEFAULT_PREDICTOR_FILES.items():
        parser.add_argument(
            f"--{key.lower().replace('_', '-')}-file",
            type=Path,
            default=default_path,
            help=f"Path to {key} input file.",
        )

    return parser.parse_args()


def main() -> None:
    """Run the final predictor-merge workflow."""
    args = parse_args()

    predictor_files_raw = {
        key: getattr(args, f"{key.lower()}_file")
        for key in DEFAULT_PREDICTOR_FILES
    }
    predictors_dir = args.predictors_dir
    predictor_files = {
        key: (
            path if path.is_absolute() or predictors_dir is None else predictors_dir / path
        )
        for key, path in predictor_files_raw.items()
    }

    missing_files = [str(p) for p in [args.experimental_master_file, *predictor_files.values()] if not p.exists()]
    if missing_files:
        raise FileNotFoundError(
            "Missing required input files:\n- " + "\n- ".join(missing_files)
        )

    out_file = args.out_file
    if out_file is None:
        stem = args.experimental_master_file.stem
        out_file = args.experimental_master_file.with_name(f"{stem}{DEFAULT_OUT_FILE_SUFFIX}")

    final_df = load_experimental_master(args.experimental_master_file)
    predictors = prepare_all_predictors(predictor_files)
    final_df = merge_all_predictors(final_df, predictors)

    print_final_summary(final_df)

    out_file.parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(out_file, sep="\t", index=False)

    print("\nSaved:")
    print(out_file)


if __name__ == "__main__":
    main()
