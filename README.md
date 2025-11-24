# fastMETA
fastMETA: An ultra fast and efficient software tool for multivariate meta-analysis of GWAS

# Web tool 
Web tool available at: https://compgen.dib.uth.gr/fastMETA/

# Installation guide
fastMETA is written in Python (ver. 3.8.2)

1)	Clone or download fastMETA from: https://github.com/pbagos/fastMETA
  ```
  git clone  https://github.com/pbagos/fastMETA
  ```

2)	After downloading the .zip folder of fastMETA from GitHub, extract it to a working directory. 

3)	Το install the requirements, pip needs to be installed. Download the script for pip, from: https://bootstrap.pypa.io/get-pip.py.

4)	Open a terminal/command prompt, cd to the folder containing the get-pip.py file and run:
    ```
    python3 get-pip.py
    ```

5)	To install the mentioned requirements with pip, open a terminal/command prompt and run:
    ```
    pip install -r  requirements.txt
    ```

# Usage 
## Input file structure 

Ensure your data follow this structure. You may include more than two pairs of BETA and SE columns. The variable column can contain any identifiers—not necessarily rsIDs.

| variable | BETA1          | BETA2          | SE1             | SE2             |
|----------|----------------|----------------|------------------|------------------|
| var1     | 0.0715430377   | 0.0218272071   | 0.0220505147     | 0.0220699109     |
| var1     | 0.0966807914   | 0.0854040595   | 0.0604912844     | 0.0607718594     |
| var1     | 0.1383389831   | 0.0450079846   | 0.0998447463     | 0.1000947719     |
| var2     | 0.1298049745   | 0.1219084187   | 0.0261241204     | 0.0261034413     |
| var2     | 0.0737540951   | 0.0285241921   | 0.0219609195     | 0.0219791685     |
| var3     | 0.0912965958   | 0.0744658800   | 0.0605097372     | 0.0607870714     |
| var3     | 0.1513790485   | 0.0512020952   | 0.0978128888     | 0.0980551991     |

## Execution command 
```
python3 fastmeta.py --method [method] --input_file [input_file_name] --output_file [output_file_name]  --het_est [heterogeneity_estimator]
```
## Arguments 
 
`fastMETA` accepts the following command-line arguments:

- `--input_file`: Input tab-separated file with columns: `variable`, `BETA1`, `SE1`, `BETA2`, `SE2`, ...
- `--output_file`: Output tab-separated file to save the meta-analysis results.
- `--het_est`: Heterogeneity estimator to use. Options:
  - `'DL'` (DerSimonian–Laird)
  - `'ANOVA'` (Cochran–ANOVA)
  - `'SJ'` (Sidik–Jonkman)
  - `'FE'` (Fixed Effects Only)  
  **Note:** `--het_est` Only applicable for `method2` and `method3`.
- `--method`: Select the multivariate meta-analysis method. Options:
  - `'method1'` — xmeta MMoM method (no `het_est` required)
  - `'method2'` — computes correlation between beta coefficients of studies for each variable
  - `'method3'` — assumes one shared correlation coefficient across all variables

### Additional arguments for `method3`

The following options are **only used when** `--method method3` is selected and are ignored for `method1` and `method2`.

- `--z_cut_off FLOAT`  
  Absolute Z-score threshold used to filter rows when constructing the correlation matrix.  
  Only rows where **all** non-missing |Z| values across studies are ≤ `z_cut_off` are retained for the correlation matrix.  
  If not provided, no Z-based filtering is applied.  
  *Default:* no Z-based filtering.

- `--R2_cut_off`  
  Enable LD-based pruning using TOP-LD for the correlation matrix. Download TOP-LD reference panel in parquet files [here](http://195.251.108.185/ref_panels/TOP_LD/EUR/SNV/)  
  When this flag is set, TOP-LD is used across chromosomes to identify `variable` entries that are in pairs with LD  
  `R² >= ld_R2_threshold`, and these variants are excluded from the correlation matrix.  
  This affects **only** the correlation matrix estimation, not the meta-analysis step itself.  
  *Default:* flag not set (no LD-pruning is performed).

- `--ld_population CODE`  
  Population code for TOP-LD paths ( `EUR`, `AFR`, `EAS` and `SAS`).  
  This is **required** for `method3` if `--R2_cut_off` is set, so that the appropriate LD reference panel is used.  
  *Default:* not set.

- `--ld_maf_threshold FLOAT`  
  Minor allele frequency (MAF) threshold applied when reading the TOP-LD MAF file in `method3`.  
  Variants with MAF below this threshold are excluded from the LD-based filtering step.  
  *Default:* `0.0`.

- `--ld_R2_threshold FLOAT`  
  LD `R²` threshold used when applying TOP-LD–based filtering in `method3`.  
  Variant pairs with `R² >= ld_R2_threshold` are considered in LD, and their corresponding `variable` entries are excluded  
  from the correlation matrix if `--R2_cut_off` is enabled.  
  *Default:* `0.2`.

## Example
```bash
python3 fastmeta.py --method method3 --input_file example_input.txt --output_file results.txt --het_est DL
```
 ## Example
```
 python3 fastmeta.py --method method3  --input_file example_input.txt  --output_file results.txt --het_est DL
```
