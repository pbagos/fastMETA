# fastMETA
fastMETA: an ultra fast and efficient software tool for multivariate meta-analysis of GWAS

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
    python get-pip.py
    ```

5)	To install the mentioned requirements with pip, open a terminal/command prompt and run:
    ```
    pip install -r  requirements.txt
    ```
    
# Execution command 
```
python fastmeta.py --method [method] --input_file [input_file_name] --output_file [output_file_name]  --het_est [heterogeneity_estimator]
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
  **Note:** Only applicable for `method1` and `method2`.
- `--method`: Select the multivariate meta-analysis method. Options:
  - `'method1'` — xmeta MMoM method (no `het_est` required)
  - `'method2'` — computes correlation between beta coefficients of studies for each variable
  - `'method3'` — assumes one shared correlation coefficient across all variables

 
