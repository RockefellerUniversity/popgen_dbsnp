[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14367749.svg)](https://doi.org/10.5281/zenodo.14367749)

## PopGen dbSNP

use the dbSNP bzip2 json file and a gtf as input, calculate the Tajima's D per gene over the exon region.


## Usage:

<code>
for num in $(seq 1 22) ; 
do 
./target/release/popgen_dbSNP --gtf ../chr${num}.gtf --dbsnp ../refsnp-chr${num}.json.bz2 --bin 0 --studyname 1000Genomes > chr${num}.popgen_1kG.tsv 2>chr${num}.popgen_1kG.err & ; 
done
</code>
