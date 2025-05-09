#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import sys
import argparse
from scipy.stats import chi2

def mmom_multi(ys, vars_):
 
    # Extract effect sizes
    ys1 = ys[:, 0]
    ys2 = ys[:, 1]

    # Extract variances
    vars1 = vars_[:, 0]
    vars2 = vars_[:, 1]

    # Inverse-variance weights
    w1 = 1.0 / vars1
    w2 = 1.0 / vars2

    # Weighted means
    y1_weight = np.sum(w1 * ys1) / np.sum(w1)
    y2_weight = np.sum(w2 * ys2) / np.sum(w2)

    # Count number of studies not imputed (i.e., not missing)
    n1 = np.sum(vars1 < 1e4)
    n2 = np.sum(vars2 < 1e4)

    # Q statistics
    Q1 = np.sum(w1 * (ys1 - y1_weight) ** 2)
    Q2 = np.sum(w2 * (ys2 - y2_weight) ** 2)

    
        
    # Between-study variances (tau^2 estimates)
    tau1_2_hat = max(0, (Q1 - (n1 - 1)) / (np.sum(w1) - np.sum(w1**2) / np.sum(w1)))
    tau2_2_hat = max(0, (Q2 - (n2 - 1)) / (np.sum(w2) - np.sum(w2**2) / np.sum(w2)))

    # Adjusted weights
    w1_star = 1.0 / (vars1 + tau1_2_hat)
    w2_star = 1.0 / (vars2 + tau2_2_hat)

    # Final beta estimates
    beta1_hat = np.sum(w1_star * ys1) / np.sum(w1_star)
    beta2_hat = np.sum(w2_star * ys2) / np.sum(w2_star)

    # Variance of beta estimates
    var_beta1_hat = 1.0 / np.sum(w1_star)
    var_beta2_hat = 1.0 / np.sum(w2_star)

    # Covariance between beta1 and beta2
    w1_norm = w1_star / np.sum(w1_star)
    w2_norm = w2_star / np.sum(w2_star)
    mycov_beta = np.sum(w1_norm * w2_norm * (ys1 - beta1_hat) * (ys2 - beta2_hat))

    # Combine results
    beta_hat = np.array([beta1_hat, beta2_hat])
    sigma_hat = np.array([
        [var_beta1_hat, mycov_beta],
        [mycov_beta, var_beta2_hat]
    ])

    return {
        'beta_hat': beta_hat,
        'beta_cov': sigma_hat
    }



def main(input_file, output_file):
    import time 
    start= time.time()
    # Read the input tab-separated file.
    df = pd.read_csv(input_file, sep="\t")
    
    # Order by the 'variable' column.
    df = df.sort_values(by='variable')
    
    # Identify BETA and SE columns using regex (e.g., BETA1, BETA2, ... and SE1, SE2, ...).
    beta_cols = sorted([col for col in df.columns if re.match(r'BETA\d+', col)],
                       key=lambda x: int(re.findall(r'\d+', x)[0]))
    se_cols = sorted([col for col in df.columns if re.match(r'SE\d+', col)],
                     key=lambda x: int(re.findall(r'\d+', x)[0]))
    
    # Check that each BETA column has a matching SE column.
    if len(beta_cols) != len(se_cols):
        print("Error: The number of BETA columns does not match the number of SE columns.")
        sys.exit(1)
    
    results = []
    # Group the dataframe by 'variable'.
    grouped = df.groupby('variable')
    
    for var, group in grouped:
        row = {'variable': var}
        
        # Build matrices for effect estimates (ys) and standard errors.
        # Each column corresponds to one method.
        ys_mat = np.column_stack([group[col].values for col in beta_cols])
        se_mat = np.column_stack([group[col].values for col in se_cols])
        # Convert standard errors to variances.
        vars_mat = se_mat ** 2
        
        # Run the multivariate meta-analysis (MMoM_multi equivalent).
        res = mmom_multi(ys_mat, vars_mat)
        beta_hat = res['beta_hat']
        beta_cov = res['beta_cov']
        #print(beta_hat)
       # print(beta_cov)
        # Save results for each method: meta_BETA and meta_SE (sqrt of variance).
        for k, (bh, var_bh) in enumerate(zip(beta_hat, np.diag(beta_cov)), start=1):
            row[f'meta_BETA{k}'] = bh
            row[f'meta_SE{k}'] = np.sqrt(var_bh)
        
        # Compute the Wald statistic and associated p-value.
        try:
            inv_cov = np.linalg.pinv(beta_cov)
            # Wald statistic: beta_hat' * inv_cov * beta_hat.
            wald = beta_hat.T @ inv_cov @ beta_hat
            p_value = chi2.sf(wald, df=beta_cov.shape[0])
        except np.linalg.LinAlgError:
            wald = np.nan
            p_value = np.nan
        
        row["Wald"] = wald
        row["P"] = p_value
        
        results.append(row)
    
    # Create a DataFrame for the results and save it to a tab-separated file.
    out_df = pd.DataFrame(results)
    end = time.time()
    print(f"Executed in : {end-start} seconds.")

    out_df.to_csv(output_file, index=False, sep="\t")
    print(f"Meta-analysis results have been saved to: {output_file}")


if __name__ == '__main__':
    version = '1.0.0'
    # Header for the program.
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
    print("FAST-META: Fast Multivariate Meta-Analysis (MMoM_multi version)")
    print("Version " + version + "; April 2025")
    print("Copyright (C) 2025 Pantelis Bagos")
    print("Freely distributed under the GNU General Public Licence (GPLv3)")
    print("---------------------------------------------------------------------------------")

    parser = argparse.ArgumentParser(
        description='Multivariate random-effects meta-analysis for beta and SE pairs in a tab-separated file')
    parser.add_argument('--input_file',
                        help='Input tab-separated file with columns: variable, BETA1, SE1, BETA2, SE2, ...')
    parser.add_argument('--output_file', help='Output tab-separated file to save the meta-analysis results')
    args = parser.parse_args()
    main(args.input_file, args.output_file)
 
    