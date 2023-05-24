# Introduction
Bash scripts used in GAEP table 1.

# Running
```bash
#GAEP
gaep bkp -r hg38_simu.fasta -i hg38_simu.fastq -x pb -t 10

#Inspector
pyhton inspector.py -c hg38_simu.fasta -r hg38_simu.fastq --datatype clr -t 10 -o inspector/test_out/ 

#QUAST
pyhton quast.py -r hg38_ref.fasta -o quast_output -t 10 hg38_simu.fasta

```