######### Guideline for RAxML-NG ML tree reconstruction


### 1. Run jModelTest to estimate the nt substitution model

### 2. Run trimal to clean the alignment - no need for a job

trimal -in align_prodigals.fasta -out trimal.fasta -automated1

#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -A lp_jm_virome_group
#PBS -o stdout.$line
#PBS -e stderr.$line

module load Java
cd $VSC_SCRATCH/PhD/Run1/RVAB_trees/vp1
java -jar $VSC_DATA/softwares/jmodeltest-2.1.10/jModelTest.jar -d VP1.fas -s 11 -i -g 4 -f -BIC -a -tr 36


### 2. Run the sanity check to assess the alignment (can check the duplicates but in the end I did not do it)


module load RAxML-NG
raxml-ng-mpi --check --msa VP1.fas --model TPM1uf+I+G --prefix vp4


#Output is V1.raxml.phy

### 3. Optional: Check the number of bootstrapping for convergence (but if you go with 1000 or autoMRE, you do not need to do this step)
### Max no. of threads 5 or so. Not the more the merrier, if you specify more it gives an error. Default seed number is the current time so no need to specify actually. 


#raxml-ng-mpi --bootstrap --msa V1.raxml.phy --model TIM1+I+G --prefix V2 --seed 333 --threads 2 --bs-trees 200
#raxml-ng-mpi --bsconverge --bs-trees V2.raxml.bootstraps --prefix V3 --seed 2 --threads 2 --bs-cutoff 0.02


### 4. Run the ML tree

#raxml-ng-mpi --search --msa V1.raxml.phy --model TIM1+I+G --prefix V2 --threads 2


### OR do 3 and 4 together if you are sure about your bootstrap number or use bootstopping so that 
### RAxML decides itself how many bootstrap replicates should be run

raxml-ng-mpi --all --msa V1.raxml.phy --model TrN+I+G --prefix V2 --threads 2 --bs-trees autoMRE (can add --tree pars{25},rand{25} or --redo)


### 5. Computing branch support (no need for a job, already done in -all option!)

#raxml-ng-mpi --support --tree vp1.raxml.bestTree --bs-trees n --prefix vp1_1

### We can also check if the resulting trees differ topologically, Number of unique topologies in this tree set should be 1
### Not really the case with RVs, I have 20 instead.

raxml-ng-mpi --rfdist --tree vp1.raxml.mlTrees --prefix RF


### output is v1.raxml.support that you can visualize in FigTree or any other tree visualising software.
### I had to reroot the tree to R1 clade in order to obtain an OK looking tree. Also decreased the scale bar manually and chose bootstrap values
### to be displayed manually. Also ordered the nodes in an increasing order for better visualisation. 
### Expansion a bit, root length a bit. 
### Illustrator: Fonts: Arial. Nodes 11pt, bootstraps 7pt. Arrange the bootstrap value positions and delete the ones less than 70. 


#Guideline for RAxML-NG

#Run jModelTest to estimate the nt substitution model



#Run the sanity check to assess the alignment

module load RAxML-NG
raxml-ng-mpi --check --msa VP1.fas --model output_jmodeltest --prefix V1






########OR TO RUN JUST ONE COMMAND FINALLY

#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l walltime=12:00:00
#PBS -A lp_jm_virome_group
#PBS -m abe
#PBS -M ceren.simsek@kuleuven.be

cd $VSC_SCRATCH/PhD
module load RAxML-NG
raxml-ng-mpi --all --msa adeno_align_mafft_yeni.fas --model GTR+I+G --prefix adeno2 --threads 2 --bs-trees autoMRE
