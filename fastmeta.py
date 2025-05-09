import argparse
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(
        description='Random or Fixed Effects Meta-analysis for beta and SE pairs in a tab-separated file')
    
    parser.add_argument('--method', required=True,
                        choices=['method1', 'method2', 'method3'],
                        help="Select the multivariate meta-analysis method:\n"
                             "'method1' - xmeta MmoM method (no het_est required)\n"
                             "'method2' - correlation between betas of studies for each variable\n"
                             "'method3' - one correlation coefficient for all variables")
    
    parser.add_argument('--input_file', required=True,
                        help='Input tab-separated file with columns: variable, BETA1, SE1, BETA2, SE2, ...')
    
    parser.add_argument('--output_file', required=True,
                        help='Output tab-separated file to save the meta-analysis results')
    
    parser.add_argument('--het_est', default="DL", choices=["DL", "ANOVA", "SJ", "FE"],
                        help="Heterogeneity estimator: 'DL' (DerSimonian-Laird), "
                             "'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman), or 'FE' (Fixed Effects Only).\n"
                             "NOTE: Only applicable for method1 and method2.")

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

    # Debug: print the full command
    print(f"Running command: {' '.join(command)}")

    result = subprocess.run(command)
    
    if result.returncode != 0:
        print(f"Error: {script_to_run} failed with return code {result.returncode}")
        sys.exit(result.returncode)

if __name__ == '__main__':
    main()
