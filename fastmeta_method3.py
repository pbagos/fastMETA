#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import sys
import argparse
from scipy.stats import chi2

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
    Q = np.sum(weights_fixed * (betas - fixed_effect)**2)
    k = len(betas)
    df = k - 1
    denom = np.sum(weights_fixed) - np.sum(weights_fixed**2) / np.sum(weights_fixed)

    # Heterogeneity estimator selection
    if heterogeneity.upper() == "DL":
        tau2 = max(0, (Q - df) / denom) if denom > 0 else 0.0

    elif heterogeneity.upper() == "ANOVA":
        numerator = np.sum((betas - np.mean(betas)) ** 2) / (k - 1)
        within_variance_avg = np.sum(variances) / k
        tau2 = max(0, numerator - within_variance_avg)

    elif heterogeneity.upper() == "SJ":
        
        yi = betas
        vi = ses**2              

        # Number of studies 
        k = yi.size

      
        X = np.ones((k, 1))
        p = X.shape[1]            

        
        Y_bar  = yi.mean()        # sample mean of betas
        ymci   = yi - Y_bar       # mean centered effects
       
        tau2_0 = np.var(ymci, ddof=0)

     
      
        wi = 1.0 / (vi + tau2_0)
        W  = np.diag(wi)          # (k × k)

        
        XtWX     = X.T @ W @ X
        inv_XtWX = np.linalg.inv(XtWX)

       
        P = W - W @ X @ inv_XtWX @ X.T @ W

        beta_hat = inv_XtWX @ (X.T @ W @ yi)   # shape (1,)
        
        fitted = X @ beta_hat                  # shape (k,)
        Ymc    = yi - fitted                   # residuals
        
        RSS = float(Ymc.T @ P @ Ymc)

        tau2 = tau2_0 * RSS / (k - p)

    elif heterogeneity.upper() == "FE":
        tau2 = 0.0

    else:
        raise ValueError(f"Unsupported heterogeneity estimator: {heterogeneity}. Use 'DL', 'ANOVA', 'SJ', or 'FE'.")

    weights_random = 1.0 / (variances + tau2)
    overall_effect = np.sum(weights_random * betas) / np.sum(weights_random)
    overall_se = np.sqrt(1.0 / np.sum(weights_random))

    return overall_effect, overall_se


def main(input_file, output_file, heterogeneity):
    df = pd.read_csv(input_file, sep="\t", encoding='ISO-8859-1')
    df = df.sort_values(by='variable')

    beta_cols = sorted([col for col in df.columns if re.match(r'BETA\d+', col)],
                       key=lambda x: int(re.findall(r'\d+', x)[0]))
    se_cols = sorted([col for col in df.columns if re.match(r'SE\d+', col)],
                     key=lambda x: int(re.findall(r'\d+', x)[0]))

    R = df[beta_cols].corr(method='pearson')
    R = np.array(R)

    if len(beta_cols) != len(se_cols):
        print("Error: The number of BETA columns does not match the number of SE columns.")
        sys.exit(1)

    results = []
    grouped = df.groupby('variable')

    for var, group in grouped:
        row = {'variable': var}
        for beta_col, se_col in zip(beta_cols, se_cols):
            betas = group[beta_col].tolist()
            ses = group[se_col].tolist()
            meta_effect, meta_se = random_effects_meta_analysis(betas, ses, heterogeneity)
            row[f'meta_{beta_col}'] = meta_effect
            row[f'meta_{se_col}'] = meta_se
        results.append(row)

    out_df = pd.DataFrame(results)

    beta_cols = [col for col in out_df.columns if col.startswith("meta_BETA")]
    se_cols = [col for col in out_df.columns if col.startswith("meta_SE")]

    wald_list = []
    p_list = []

    for i, row in out_df.iterrows():
        beta_array = np.array(row[beta_cols])
        se_array = np.array(row[se_cols])
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
    out_df["R"] = R[0][1]

    out_df.to_csv(output_file, index=False, sep="\t")
    print(f"Meta-analysis results have been saved to: {output_file}")


if __name__ == '__main__':
    version = '1.0.2'
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
        description='Random or Fixed Effects Meta-analysis for beta and SE pairs in a tab-separated file')
    
    parser.add_argument('--input_file', required=True,
                        help='Input tab-separated file with columns: variable, BETA1, SE1, BETA2, SE2, ...')
    
    parser.add_argument('--output_file', required=True,
                        help='Output tab-separated file to save the meta-analysis results')
    
    parser.add_argument('--het_est', default="DL", choices=["DL", "ANOVA", "SJ", "FE"],
                        help="Heterogeneity estimator: 'DL' (DerSimonian-Laird), 'ANOVA' (Cochran-ANOVA), 'SJ' (Sidik-Jonkman), or 'FE' (Fixed Effects Only)")

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.het_est)
