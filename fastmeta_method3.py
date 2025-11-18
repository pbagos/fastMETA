import pandas as pd
import numpy as np
import re
import sys
import argparse
import gc
from scipy.stats import chi2

# Only import polars if we actually need it (when --R2_cut_off / TOP-LD is used)
try:
    import polars as pl
except ImportError:
    pl = None


def random_effects_meta_analysis(betas, ses, het_est="DL"):
    """
    Perform a (random or fixed) effects meta-analysis using a specified heterogeneity estimator.

    Parameters:
        betas (list or array): Effect estimates from individual studies.
        ses (list or array): Standard errors corresponding to the effect estimates.
        het_est (str): Heterogeneity estimator ('DL', 'ANOVA', 'SJ', or 'FE').

    Returns:
        tuple: (overall_effect, overall_se)
    """
    betas = np.array(betas)
    ses = np.array(ses)
    mask = ~np.isnan(betas) & ~np.isnan(ses)
    betas = betas[mask]
    ses = ses[mask]
    heterogeneity = het_est

    if len(betas) == 0:
        return np.nan, np.nan

    variances = ses ** 2
    weights_fixed = 1.0 / variances
    fixed_effect = np.sum(weights_fixed * betas) / np.sum(weights_fixed)
    Q = np.sum(weights_fixed * (betas - fixed_effect) ** 2)
    k = len(betas)
    df = k - 1
    denom = np.sum(weights_fixed) - np.sum(weights_fixed ** 2) / np.sum(weights_fixed)

    # Heterogeneity estimator selection
    if heterogeneity.upper() == "DL":
        tau2 = max(0, (Q - df) / denom) if denom > 0 else 0.0

    elif heterogeneity.upper() == "ANOVA":
        numerator = np.sum((betas - np.mean(betas)) ** 2) / (k - 1)
        within_variance_avg = np.sum(variances) / k
        tau2 = max(0, numerator - within_variance_avg)

    elif heterogeneity.upper() == "SJ":
        t0 = max(0.01, np.sum((betas - np.mean(betas)) ** 2) / k)
        vi = (variances / t0) + 1
        weights_random = 1.0 / (vi + t0)
        mmwo = np.sum(weights_random * betas) / np.sum(weights_random)
        tau2 = max(0, (1 / (k - 1)) * np.sum(((betas - mmwo) ** 2) / vi))

    elif heterogeneity.upper() == "FE":
        tau2 = 0.0

    else:
        raise ValueError(
            f"Unsupported heterogeneity estimator: {heterogeneity}. "
            "Use 'DL', 'ANOVA', 'SJ', or 'FE'."
        )

    weights_random = 1.0 / (variances + tau2)
    overall_effect = np.sum(weights_random * betas) / np.sum(weights_random)
    overall_se = np.sqrt(1.0 / np.sum(weights_random))

    return overall_effect, overall_se


# -------------------------------------------------------------------------
# TOP-LD helper functions (adapted from your provided code)
# -------------------------------------------------------------------------
def TOP_LD_info(rs_list, chrom, population, maf_threshold, R2_threshold, imp_snp_list=None):
    """
    Query TOP-LD Parquet files and return LD pairs with R2 >= R2_threshold.

    Returns a pandas DataFrame with columns:
        pos1, pos2, R2, +/-corr, Dprime, rsID1, rsID2,
        MAF1, MAF2, REF1, ALT1, REF2, ALT2
    """
    if pl is None:
        print(
            "Error: polars is required for TOP-LD lookups but is not installed.\n"
            "Please install it with: pip install polars"
        )
        sys.exit(1)

    print(f"Loading TOP-LD files for chromosome {chrom}...")

    # 1) Combine rsIDs up-front
    if imp_snp_list:
        all_rsids = list(set(rs_list) | set(imp_snp_list))
    else:
        all_rsids = rs_list

    # 2) Lazy-scan the MAF file, project only needed cols, then filter
    maf_path = (
        f"ref_panels/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_info_annotation.parquet"
    )
    maf_lazy = (
        pl.scan_parquet(maf_path)
        .select(["Position", "rsID", "MAF", "REF", "ALT"])
        .filter(pl.col("MAF") >= maf_threshold)
        .filter(pl.col("rsID").is_in(all_rsids))
    )

    # 3) Lazy-scan the LD file, project and filter by R²
    ld_path = (
        f"ref_panels/TOP_LD/{population}/SNV/"
        f"{population}_chr{chrom}_no_filter_0.2_1000000_LD.parquet"
    )
    ld_lazy = (
        pl.scan_parquet(ld_path)
        .select(["SNP1", "SNP2", "R2", "+/-corr", "Dprime"])
        .filter(pl.col("R2") >= R2_threshold)
    )

    # 4) Prepare two MAF views for joining
    maf1 = maf_lazy.rename({
        "Position": "SNP1", "rsID": "rsID1", "MAF": "MAF1",
        "REF": "REF1", "ALT": "ALT1"
    })
    maf2 = maf_lazy.rename({
        "Position": "SNP2", "rsID": "rsID2", "MAF": "MAF2",
        "REF": "REF2", "ALT": "ALT2"
    })

    # 5) Join LD ↔ MAF1 ↔ MAF2 all lazily
    joined = (
        ld_lazy
        .join(maf1, on="SNP1", how="inner")
        .join(maf2, on="SNP2", how="inner")
    )

    # 6) If you provided an imp_snp_list, filter rsID roles
    if imp_snp_list:
        joined = joined.filter(
            (pl.col("rsID1").is_in(rs_list)) &
            (pl.col("rsID2").is_in(imp_snp_list))
        )

    # 7) Select + rename final output columns
    final_lazy = joined.select([
        "rsID1", "rsID2"

    ])

    # 8) Execute once in streaming mode (very memory-efficient)
    result = final_lazy.collect()

    if result.is_empty():
        print(f"TOP-LD: No SNPs found above the given thresholds on chromosome {chrom}.")

    # 9) Cleanup
    del maf_lazy, maf1, maf2, ld_lazy, joined, final_lazy
    gc.collect()

    return result.to_pandas()


def TOP_LD_process(study_df, r2threshold, population, maf_input, chromosome, imp_snp_list=None):
    """
    Wrapper that takes a DataFrame with a 'SNP' column and
    returns LD pairs from TOP-LD for those SNPs for one chromosome.
    """
    outputData = TOP_LD_info(
        list(study_df['SNP']),
        chromosome,
        population,
        maf_input,
        r2threshold,
        imp_snp_list
    )
    return outputData


# -------------------------------------------------------------------------
# MAIN META-ANALYSIS FUNCTION
# -------------------------------------------------------------------------
def main(input_file,
         output_file,
         heterogeneity,
         z_cut_off=None,
         r2_cut_off=False,
         r2_parquet="unique_rsIDs_chr1_22.parquet",  # kept for backward compatibility, not used now
         ld_population=None,
         ld_chr=None,
         ld_maf_threshold=0.0,
         ld_R2_threshold=0.2):

    df = pd.read_csv(input_file, sep="\t", encoding='ISO-8859-1')
    df = df.sort_values(by='variable')

    beta_cols = sorted(
        [col for col in df.columns if re.match(r'BETA\d+', col)],
        key=lambda x: int(re.findall(r'\d+', x)[0])
    )
    se_cols = sorted(
        [col for col in df.columns if re.match(r'SE\d+', col)],
        key=lambda x: int(re.findall(r'\d+', x)[0])
    )

    if len(beta_cols) != len(se_cols):
        print("Error: The number of BETA columns does not match the number of SE columns.")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # OPTIONAL Z-BASED FILTERING FOR CORRELATION MATRIX
    # -------------------------------------------------------------------------
    if z_cut_off is not None:
        with np.errstate(divide='ignore', invalid='ignore'):
            beta_matrix = df[beta_cols].to_numpy(dtype=float)
            se_matrix = df[se_cols].to_numpy(dtype=float)
            z_matrix = beta_matrix / se_matrix  # element-wise

        # Valid if |Z| <= z_cut_off OR Z is NaN (we ignore missing values)
        valid_z = np.logical_or(
            np.isnan(z_matrix),
            (z_matrix >= -z_cut_off) & (z_matrix <= z_cut_off)
        )

        # We keep rows where ALL non-missing Zs are within [-z_cut_off, z_cut_off]
        valid_rows_for_corr = np.all(valid_z, axis=1)
        df_for_corr = df.loc[valid_rows_for_corr].copy()
    else:
        # No Z-based filtering: use all rows for the correlation matrix
        df_for_corr = df.copy()

    # -------------------------------------------------------------------------
    # OPTIONAL R2-BASED EXCLUSION USING TOP-LD ACROSS ALL CHROMOSOMES
    #
    # New behaviour:
    #   If r2_cut_off is True:
    #     - Build rsID list from BOTH 'variable' and 'rsID' columns
    #       (variable only used when it looks like an rsID, i.e. '^rs\\d+$').
    #     - Loop over chromosomes:
    #          * If ld_chr is None or 'all' -> 1..22
    #          * Else ld_chr can be a comma-separated list, e.g. '1,2,10'
    #     - Query TOP-LD per chromosome, collect all LD pairs with R2 >= ld_R2_threshold.
    #     - Exclude from df_for_corr any row whose 'variable' OR 'rsID'
    #       is in the union of rsID1 across chromosomes.
    #
    # The meta-analysis (df) is NOT altered.
    # -------------------------------------------------------------------------
    if r2_cut_off:
        if pl is None:
            print(
                "Error: polars is required for --R2_cut_off but is not installed.\n"
                "Please install it with: pip install polars"
            )
            sys.exit(1)

        if 'variable' not in df_for_corr.columns:
            print("Error: 'variable' column is required in the input file.")
            sys.exit(1)

        if ld_population is None:
            print(
                "Error: --ld_population must be provided "
                "when --R2_cut_off is set (TOP-LD lookup)."
            )
            sys.exit(1)

        # Decide which chromosomes to use
        if (ld_chr is None) or (str(ld_chr).lower() == "all"):
            chrom_list = [str(c) for c in range(1, 23)]  # 1..22
        else:
            # allow comma-separated list, e.g. "1,2,10"
            chrom_list = [c.strip() for c in str(ld_chr).split(",") if c.strip()]

        # ---- Build rsID list from variable and rsID columns ----
        rsid_set = set()

        # # Use rsID column if present
        # if 'rsID' in df_for_corr.columns:
        #     rsid_set.update(
        #         df_for_corr['rsID']
        #         .dropna()
        #         .astype(str)
        #         .tolist()
        #     )

        # Also use any 'variable' entries that look like rsIDs (rs12345, etc.)
        var_series = df_for_corr['variable'].dropna().astype(str)


        rs_list = pd.unique(var_series)

        print(f"TOP-LD R2 filter: using {len(rs_list)} unique rsIDs from 'variable' column.")

        if len(rs_list) == 0:
            print("TOP-LD R2 filter: no rsIDs found in 'variable'; skipping LD-based exclusion.")
        else:
            # Prepare study_df with a 'SNP' column holding all rsIDs
            study_df = pd.DataFrame({'SNP': rs_list})

            ld_frames = []

            for chrom in chrom_list:
                ld_df_chr = TOP_LD_process(
                    study_df=study_df,
                    r2threshold=ld_R2_threshold,
                    population=ld_population,
                    maf_input=ld_maf_threshold,
                    chromosome=chrom,
                    imp_snp_list=None
                )

                if ld_df_chr is not None and not ld_df_chr.empty:
                    ld_frames.append(ld_df_chr)

            if ld_frames:
                ld_df_all = pd.concat(ld_frames, ignore_index=True)
            else:
                ld_df_all = pd.DataFrame()

            if ld_df_all is not None and not ld_df_all.empty:
                # All pairs already satisfy R2 >= ld_R2_threshold
                rsids_to_exclude = set(ld_df_all['rsID1'].dropna().astype(str).unique())

                before_rows = df_for_corr.shape[0]

                # Exclude if EITHER 'variable' OR 'rsID' is in rsids_to_exclude
                exclude_mask = df_for_corr['variable'].astype(str).isin(rsids_to_exclude)



                df_for_corr = df_for_corr[~exclude_mask].copy()
                after_rows = df_for_corr.shape[0]
                n_excluded = before_rows - after_rows

                # print(
                #     f"R2 filter (TOP-LD, chromosomes {','.join(chrom_list)}): "
                #     f"excluded {n_excluded} variants from correlation because "
                #     f"they appear as rsID1 in pairs with R2 >= {ld_R2_threshold}."
                # )
            else:
                print(
                    f"R2 filter (TOP-LD): no LD pairs found with R2 >= {ld_R2_threshold}; "
                    "no variants excluded from correlation."
                )

    # >>> Reporting of variable counts <<<
    n_initial_vars = df['variable'].nunique()
    n_final_vars = df_for_corr['variable'].nunique()
    print(
        f"Variables for correlation: {n_initial_vars} initially, "
        f"{n_final_vars} after filtering (Z / R2)."
    )
    # >>> End <<<

    # -------------------------------------------------------------------------
    # CORRELATION MATRIX CALCULATION
    # -------------------------------------------------------------------------
    if df_for_corr.shape[0] < 2:
        # If too few rows remain for a sensible correlation, fall back to identity
        R = np.eye(len(beta_cols))
    else:
        R = df_for_corr[beta_cols].corr(method='pearson').to_numpy()

    # -------------------------------------------------------------------------
    # META-ANALYSIS PER VARIABLE
    # -------------------------------------------------------------------------
    results = []
    grouped = df.groupby('variable')

    for var, group in grouped:
        row = {'variable': var}

        # N = number of studies (rows) where this variable appears
        study_count = group.shape[0]

        # Meta-analysis across all rows for each BETA/SE pair
        for beta_col, se_col in zip(beta_cols, se_cols):
            betas = group[beta_col].tolist()
            ses = group[se_col].tolist()

            meta_effect, meta_se = random_effects_meta_analysis(betas, ses, heterogeneity)
            row[f'meta_{beta_col}'] = meta_effect
            row[f'meta_{se_col}'] = meta_se

        row['N'] = study_count
        results.append(row)

    out_df = pd.DataFrame(results)

    beta_cols_meta = [col for col in out_df.columns if col.startswith("meta_BETA")]
    se_cols_meta = [col for col in out_df.columns if col.startswith("meta_SE")]

    # -------------------------------------------------------------------------
    # MULTIVARIATE WALD TEST
    # -------------------------------------------------------------------------
    wald_list = []
    p_list = []

    for i, row in out_df.iterrows():
        beta_array = np.array(row[beta_cols_meta], dtype=float)
        se_array = np.array(row[se_cols_meta], dtype=float)
        S_u = np.diag(se_array)

        var_cov_matrix = S_u @ R @ S_u
        try:
            inv_cov = np.linalg.pinv(np.array(var_cov_matrix, dtype=np.float64))
            wald = beta_array.T @ inv_cov @ beta_array
            p_value = chi2.sf(wald, df=beta_array.shape[0])
        except np.linalg.LinAlgError:
            wald = np.nan
            p_value = np.nan

        wald_list.append(wald)
        p_list.append(p_value)

    out_df["Wald"] = wald_list
    out_df["P"] = p_list

    # For compatibility with your original code: store one R entry (correlation between first two studies)
    if R.shape[0] >= 2 and R.shape[1] >= 2:
        out_df["R"] = R[0][1]
    else:
        out_df["R"] = np.nan

    # Ensure N is the last column in the output
    cols = list(out_df.columns)
    if "N" in cols:
        cols.append(cols.pop(cols.index("N")))
        out_df = out_df[cols]

    out_df.to_csv(output_file, index=False, sep="\t")
    print(f"Meta-analysis results have been saved to: {output_file}")


if __name__ == '__main__':
    version = '1.0.0'
    print("""
      --------------------------------------------------------------------------------
      |   ########    ###     ######  ######## ##     ## ######## ########    ###     | 
      |   ##         ## ##   ##    ##    ##    ###   ### ##          ##      ## ##    | 
      |   ##        ##   ##  ##          ##    #### #### ##          ##     ##   ##   | 
      |   ######   ##     ##  ######     ##    ## ### ## ######      ##    ##     ##  | 
      |   ##       #########       ##    ##    ##     ## ##          ##    #########  | 
      |   ##       ##     ## ##    ##    ##    ##     ## ##          ##    ##     ##  | 
      |   ##       ##     ##  ######     ##    ##     ## ########    ##    ##     ##  | 
      --------------------------------------------------------------------------------
    """)
    print("FAST-META: Fast Multivariate Meta-Analysis")
    print("Version " + version + "; April 2025")
    print("Copyright (C) 2025 Pantelis Bagos")
    print("Freely distributed under the GNU General Public Licence (GPLv3)")
    print("---------------------------------------------------------------------------------")

    parser = argparse.ArgumentParser(
        description='Random or Fixed Effects Meta-analysis for beta and SE pairs in a tab-separated file'
    )

    parser.add_argument(
        '--input_file', required=True,
        help='Input tab-separated file with columns: variable, BETA1, SE1, BETA2, SE2, ...'
    )

    parser.add_argument(
        '--output_file', required=True,
        help='Output tab-separated file to save the meta-analysis results'
    )

    parser.add_argument(
        '--het_est', default="DL", choices=["DL", "ANOVA", "SJ", "FE"],
        help=(
            "Heterogeneity estimator: 'DL' (DerSimonian-Laird), "
            "'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman), or 'FE' (Fixed Effects Only)"
        )
    )

    parser.add_argument(
        '--z_cut_off', type=float, default=None,
        help=(
            "Absolute Z threshold used to filter rows for the correlation matrix "
            "(keep rows where all non-missing |Z| <= z_cut_off). "
            "If not provided, no Z-based filtering is applied."
        )
    )

    parser.add_argument(
        '--R2_cut_off', action='store_true',
        help=(
            "If set, use TOP-LD across chromosomes to identify rsID1 variants in pairs "
            "with LD R2 >= ld_R2_threshold and exclude those variants from the "
            "correlation matrix (matching against both 'variable' and 'rsID' columns). "
            "This only affects the correlation matrix, not the meta-analysis itself."
        )
    )

    # Kept for backward compatibility, but no longer used in the new logic.
    parser.add_argument(
        '--R2_parquet', default="unique_rsIDs_chr1_22.parquet",
        help=(
            "Deprecated: previously used to pass a Parquet file with UNIQUE rsIDs to be EXCLUDED "
            "from the correlation matrix. Ignored in the current version."
        )
    )

    # New arguments for TOP-LD usage
    parser.add_argument(
        '--ld_population', default=None,
        help="Population code for TOP-LD paths (e.g. EUR, AFR, EAS). Required if --R2_cut_off is set."
    )

    parser.add_argument(
        '--ld_chr', default=None,
        help=(
            "Chromosome(s) for TOP-LD paths, e.g. '1', '1,2,3' or 'all'. "
            "If omitted or 'all', chromosomes 1-22 are used."
        )
    )

    parser.add_argument(
        '--ld_maf_threshold', type=float, default=0.0,
        help="MAF threshold applied when reading TOP-LD MAF file (default: 0.0)."
    )

    parser.add_argument(
        '--ld_R2_threshold', type=float, default=0.2,
        help="R2 threshold for LD-based exclusion via TOP-LD (default: 0.2)."
    )

    args = parser.parse_args()

    main(
        args.input_file,
        args.output_file,
        args.het_est,
        args.z_cut_off,
        args.R2_cut_off,
        args.R2_parquet,
        args.ld_population,
        args.ld_chr,
        args.ld_maf_threshold,
        args.ld_R2_threshold
    )
