We generated 500 reads of measles K01711 and 10 reads for staphylococcus (AP018562) with 
ART software single read. We concatenated the first 10 reads followed by 500
reads into a single file data.fastq

The directory krakendb is a simple databases made of only one single sequence
(K01711)

https://www.ncbi.nlm.nih.gov/nuccore/AP018562


to build DB::

    makeblastdb  -in input.fasta -dbtype nucl -out temp
