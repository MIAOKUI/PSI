## Quit Start 
### For Shell version:

### For python version: 
1  Counting inclusion count
```bash 
./psi_python/dexseq_count.py -p no -s no -r pos flattened_gff_file.gtf  alignment_file.bam basename.inclusion
```

2  Counting exclusion count 
```bash 
./psi_python/exclusion_count.py  flattened_gff_file.gtf junctions.bed basename.exclusion
```

3  Calculate the psi 
```bash 
./psi_python/psi_calculation.py flattened_gff_file.gtf basename.inclusion basename.exclusion basename
```
#### Sample Result
```shell 
ENSRNOG00000033734:001	1.0	1	44487	0	13	48872971	48873009	+	38
ENSRNOG00000033734:002	1.0	1	45933	0	13	48873964	48873975	+	11
ENSRNOG00000033734:003	0.170570528295	1	7490	32484	13	48874085	48874097	+	12
ENSRNOG00000033734:004	0.933449251441	1	51960	2843	13	48876463	48876493	+	30
```
*  Exonic part ID
*  PSI value
*  Low inclusion filter tag, 1 pass, 0.2 failed 
*  Inclusion raw count from dexseq_count.py
*  Exclusion raw count from exlcusion_count.py 
*  Chromosome ID
*  Start coordinate of exonic part 
*  End coordinate of exonic part 
*  Strand tag
*  Length of exonic part

If the set **--detail** **yes**, then only showing first three fields like following.
```shell
ENSRNOG00000033734:001	1.0	1
ENSRNOG00000033734:002	1.0	1
ENSRNOG00000033734:003	0.170570528295	1
ENSRNOG00000033734:004	0.933449251441	1
``` 
