# fastMETA
fastMETA: an ultra fast and efficient software tool for multivariate meta-analysis of GWAS

## Installation guide
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
