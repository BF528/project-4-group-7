zless -S gencode.v33.annotation.gtf.gz | grep -v "#" | awk '$3=="transcript"' | cut -f9 | tr -s ";" " " | awk '{print$4"\t"$2}' | sort | uniq |  sed 's/\"//g' >  tx2gene_ensemble.tsv




