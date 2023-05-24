# Introduction
These scripts are used to simulate the misassemblies using a template genome.

# Running
```bash
#Randomly generate the positions of misassemlies
perl simu_misassembly_posi.pl template.fasta > position.txt

#Remove the redundant positions 
sort -k1V -k2n position.txt | perl move_redundant.pl > position_redun.txt

#Introduce misassemblies by positions. Two FASTA files will be output: one for the reference and one for simulation.
perl simu_misassembly.pl template.fasta
```