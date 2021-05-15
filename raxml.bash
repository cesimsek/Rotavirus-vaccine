######### RAxML-NG ML tree reconstruction


### 1. Run jModelTest to estimate the nt substitution model

# Submit job

module load Java
cd $VSC_SCRATCH/etc
java -jar jModelTest.jar -d VP1.fas -s 11 -i -g 4 -f -BIC -a -tr 36


### 2. Run trimal to clean the alignment

trimal -in align_prodigals.fasta -out trimal.fasta -automated1


### 3. Run the sanity check to assess the alignment 

module load RAxML-NG
raxml-ng-mpi --check --msa VP1.fas --model TrN+I+G --prefix V1

### 3. Run the ML tree

#raxml-ng-mpi --search --msa V1.raxml.phy --model TrN+I+G --prefix V2 --threads 2


### OR do 3 and 4 together if you are sure about your bootstrap number or use bootstopping so that 
### RAxML decides itself how many bootstrap replicates should be run

raxml-ng-mpi --all --msa V1.raxml.phy --model TrN+I+G --prefix V1 --threads 2 --bs-trees autoMRE (can add --tree pars{25},rand{25} or --redo
