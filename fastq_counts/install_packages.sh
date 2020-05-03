conda install -c conda-forge -c bioconda \
	python=3.7.6 star=2.7.3a samtools=1.10 \
	salmon=1.2.0 bioconductor-alevinqc=1.2.0 \
	bwa=0.7.17

pip install -r requirements.txt

Rscript r_packages.R
conda list --export > conda_packages.txt 

