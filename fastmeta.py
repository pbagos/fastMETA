import argparse
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(
        description='Random or Fixed Effects Meta-analysis for beta and SE pairs in a tab-separated file'
    )

    parser.add_argument(
        '--method',
        required=True,
        choices=['method1', 'method2', 'method3'],
        help=(
            "Select the multivariate meta-analysis method:\n"
            "'method1' - xmeta MmoM method (no het_est required)\n"
            "'method2' - correlation between betas of studies for each variable\n"
            "'method3' - one correlation coefficient for all variables"
        )
    )

    parser.add_argument(
        '--input_file',
        required=True,
        help='Input tab-separated file with columns: variable, BETA1, SE1, BETA2, SE2, ...'
    )

    parser.add_argument(
        '--output_file',
        required=True,
        help='Output tab-separated file to save the meta-analysis results'
    )

    parser.add_argument(
        '--het_est',
        default="DL",
        choices=["DL", "ANOVA", "SJ", "FE"],
        help=(
            "Heterogeneity estimator: 'DL' (DerSimonian-Laird), "
            "'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman), or 'FE' (Fixed Effects Only).\n"
            "NOTE: Only applicable for method2 and method3 (method1/xmeta does not use it)."
        )
    )

    # -----------------------------
    # method3-specific arguments
    # -----------------------------
    parser.add_argument(
        '--z_cut_off',
        type=float,
        default=None,
        help=(
            "Absolute Z threshold used to filter rows for the correlation matrix "
            "(keep rows where all non-missing |Z| <= z_cut_off). "
            "If not provided, no Z-based filtering is applied. "
            "Only used by method3."
        )
    )

    parser.add_argument(
        '--R2_cut_off',
        action='store_true',
        help=(
            "If set (for method3), use TOP-LD across chromosomes to identify rsID1 variants in "
            "pairs with LD R2 >= ld_R2_threshold and exclude those variants from the correlation "
            "matrix (matching against 'variable'). This only affects the correlation matrix, "
            "not the meta-analysis itself."
        )
    )

    # Kept for backward compatibility with method3’s script
    parser.add_argument(
        '--R2_parquet',
        default="unique_rsIDs_chr1_22.parquet",
        help=(
            "Deprecated in the current method3 logic but kept for backward compatibility. "
            "Previously used to pass a Parquet file with UNIQUE rsIDs to be excluded from the "
            "correlation matrix."
        )
    )

    parser.add_argument(
        '--ld_population',
        default=None,
        help=(
            "Population code for TOP-LD paths (e.g. EUR, AFR, EAS). "
            "Required for method3 if --R2_cut_off is set."
        )
    )

    parser.add_argument(
        '--ld_chr',
        default=None,
        help=(
            "Chromosome(s) for TOP-LD paths used by method3, e.g. '1', '1,2,3' or 'all'. "
            "If omitted or 'all', chromosomes 1-22 are used."
        )
    )

    parser.add_argument(
        '--ld_maf_threshold',
        type=float,
        default=0.0,
        help=(
            "MAF threshold applied when reading TOP-LD MAF file in method3 (default: 0.0)."
        )
    )

    parser.add_argument(
        '--ld_R2_threshold',
        type=float,
        default=0.2,
        help=(
            "R2 threshold for LD-based exclusion via TOP-LD in method3 (default: 0.2)."
        )
    )

    args = parser.parse_args()

    method_map = {
        'method1': 'fastmeta_method1.py',
        'method2': 'fastmeta_method2.py',
        'method3': 'fastmeta_method3.py'
    }

    script_to_run = method_map.get(args.method)
    if not script_to_run:
        print(f"Unknown method: {args.method}")
        sys.exit(1)

    # Basic command
    command = [
        'python3', script_to_run,
        '--input_file', args.input_file,
        '--output_file', args.output_file
    ]

    # Add het_est only for method2 and method3
    if args.method in ['method2', 'method3']:
        command.extend(['--het_est', args.het_est])

    # Add method3-specific options
    if args.method == 'method3':
        # Z-based filtering for correlation matrix
        if args.z_cut_off is not None:
            command.extend(['--z_cut_off', str(args.z_cut_off)])

        # R2-based LD filtering
        if args.R2_cut_off:
            command.append('--R2_cut_off')

        # R2_parquet (backwards compatibility)
        if args.R2_parquet is not None:
            command.extend(['--R2_parquet', args.R2_parquet])

        # TOP-LD related parameters
        if args.ld_population is not None:
            command.extend(['--ld_population', args.ld_population])

        if args.ld_chr is not None:
            command.extend(['--ld_chr', args.ld_chr])

        if args.ld_maf_threshold is not None:
            command.extend(['--ld_maf_threshold', str(args.ld_maf_threshold)])

        if args.ld_R2_threshold is not None:
            command.extend(['--ld_R2_threshold', str(args.ld_R2_threshold)])

    # Debug: print the full command
    #print(f"Running command: {' '.join(command)}")

    result = subprocess.run(command)

    if result.returncode != 0:
        print(f"Error: {script_to_run} failed with return code {result.returncode}")
        sys.exit(result.returncode)


if __name__ == '__main__':
    main()
