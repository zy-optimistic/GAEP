# Introduction
These scripts are used to simulate the misassemblies using a template genome.

# Running
```bash
#At first, we recommend to remove the contigs less than 500k in template.fasta.

#Index
samtools faidx template.fasta

#Randomly generate the positions of misassemlies
perl simu_misassembly_posi.pl template.fasta > position.txt

#Remove the redundant positions 
sort -k1V -k2n position.txt | perl move_redundant.pl > position_redun.txt

#Introduce misassemblies by positions. Two FASTA files will be output: one for the reference and one for simulation.
perl simu_misassembly.pl template.fasta

#PacBio reads simulation
pbsim hg38.ref.fasta --prefix simu_pb --depth 50 --length-min 5000 --length-max 50000 --hmm_model XXX/P6C4.model --length-mean 20000
```
