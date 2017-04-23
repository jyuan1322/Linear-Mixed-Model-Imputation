# Linear-Mixed-Model-Imputation
Perform linear mixed model regression and soft-impute imputation in one step, minimizing error from both operations simultaneously.

Install:
pip install pyvcf
pip install pysam
sudo apt-get install tabix

Setup:
bgzip file.vcf
tabix -p vcf file.vcf.gz 

