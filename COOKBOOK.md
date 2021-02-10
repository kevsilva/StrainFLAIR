# StrainFLAIR cookbook

- - - -
**Prerequisite**

Datasets for this cookbook are available in the folder `data/minimal_example/`. It contains 4 complete genomes (fasta files) and a mixture of simulated reads with 0.001% of errors composed of 67% of reads from D4 and 33% of reads from LM33 (fastq file).

**Full commands:** For each example, the full commands are proposed at the end of the section.

**Computation time** is approximatel 3 minutes for the indexation step, 3 minutes for the mapping step and 1 minute for the query step.
- - - -

Consider 4 reference genomes:

**First**, create a *file of file* (fof) used as input for StrainFLAIR. Each line of the input fof contains exactly one fasta name. We may create this fof as follows:

```bash
ls data/minimal_example/*.fasta > list_fasta.txt
```

**Second**, run StrainFLAIR indexation step, using the input fof:

```bash
./StrainFLAIR.sh index -i list_fasta.txt -o myproject
```

The reference graph and its associated additionnal files have been created in the new directory `myproject`.

**Third**, query the variation graph by mapping any reads on it. Reads can be compressed in a fastq.gz format. 

```bash
./StrainFLAIR.sh query -g myproject/graphs/all_graphs -f1 data/minimal_example/mixture_D4_LM33_67_33.fastq -t 24 -p myproject/graphs/dict_clusters.pickle -d myproject -o minimalexample
```

That's it! Mapping files are available in `mapping/`, and abundance tables are available in `results/`.

The final result is a csv table containing each reference genome in line identified by their accession number. Columns contained the proportion of detected genes and the estimated strain-level abundance according to different computation methods. Here we found the initial ratio of around 65-68% for D4 (CP010143.1) and around 31-35% for LM33 (NZ_LN874954.1).

Here is a small example:
```
,detected_genes,mean_abund,mean_abund_nz,median_abund,median_abund_nz
NZ_LN874954.1,0.9774620284174425,31.897899712576486,34.90845134147716,31.76353773804796,34.8974830123928
CP014492.1,0.21266968325791855,0.0,0.0,0.0,0.0
CP010143.1,0.9890345649582837,68.10210028742353,65.09154865852284,68.23646226195204,65.10251698760719
AP022815.1,0.07082211638020294,0.0,0.0,0.0,0.0
```


**Full commands:**

```
ls data/minimal_example/*.fasta > list_fasta.txt
./StrainFLAIR.sh index -i list_fasta.txt -o myproject
./StrainFLAIR.sh query -g myproject/graphs/all_graphs -f1 data/minimal_example/mixture_D4_LM33_67_33.fastq -t 24 -p myproject/graphs/dict_clusters.pickle -d myproject -o minimalexample
```
