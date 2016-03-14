## Quit Start 
### For Shell version:
### For python version: 
1. Counting inclusion count
	```bash 
	./psi_python/dexseq_count.py -p no -s no -r pos flattened_gff_file.gtf  alignment_file.bam basename.inclusion
	```
2. Counting exclusion count 
	```bash 
	./psi_python/exclusion_count.py  flattened_gff_file.gtf junctions.bed basename.exclusion
	```
3. Calculate the psi 
	```bash 
	./psi_python/psi_calculation.py flattened_gff_file.gtf basename.inclusion basename.exclusion basename
	```

