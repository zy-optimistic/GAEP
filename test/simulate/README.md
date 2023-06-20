# Introduction
These scripts are used to simulate the misassemblies using a template genome.

# Running
```bash
#At first, we recommend to remove the contigs less than 500k in template.fasta and split the template.fasta from Ns.

#Index
samtools faidx template.fasta

#Randomly generate the positions of misassemlies
perl simu_misassembly_posi.pl template.fasta > position.txt

#Remove the redundant positions 
sort -k1V -k2n position.txt | perl move_redundant.pl > position_redun.txt

#Introduce misassemblies by positions. Two FASTA files will be output: one for the reference and one for simulation.
perl simu_misassembly.pl template.fasta

#PacBio reads simulation
pbsim *.ref.fasta --prefix simu_pb --depth 50 --length-min 5000 --length-max 50000 --hmm_model XXX/P6C4.model --length-mean 20000
```

# Output 
The positions of simulated misassemblies can be found in STDOUT.
```
#Format:
No  strand  contig  type  length  start  end
41	1	chr1	ins	48762	14640630	14640631
41	2	chr1	ins	48762	14697478	14746240

#For strand:
#1: Misassemblies corresponding to positions in the reference sequence.
#2: Misassembly positions in the simulated sequence.
```
Typically, retaining only the rows with a strand value of 2 will provide the positional information of the final simulated misassemblies.
