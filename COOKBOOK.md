# StrainFLAIR cookbook

- - - -
**Prerequisite**

Datasets for this cookbook are available in the folder `data/manuscript_simulated_datasets/`. It contains 8 complete genomes (fasta files) of *Escherichia coli* and several mixtures of simulated reads with 0.1% of errors as described in the manuscript (fastq.gz files).

For this cookbook, only the dataset *mix_o104_iai39_k12_100k_e0.001.fastq.gz* is used. It corresponds to a mix of 50% of *E. coli* O104:H4 2011C-3493, 33% of IAI39 and 17% of K-12 MG1655. All commands presented can be used for the other datasets, by simply changing the name of the fastq file used as shown in the third step below.

**Full commands:** The full commands are proposed at the end of the section.

**Computation time** is approximately 3 minutes for the indexation step, 30 minutes for the mapping step and 1 minute for the query step.
- - - -

Consider 7 reference genomes:

**First**, create a *file of file* (fof) used as input for StrainFLAIR. Each line of the input fof contains exactly one fasta name. The folder contains 8 genomes, however the strain BL21DE3 is only used as an unknown strain be to queried and is not indexed. We may create this fof as follows:

```bash
find data/manuscript_simulated_dataset/references/* -not -name "*bl21*" > data/list_fasta_sim.txt
```

**Second**, run StrainFLAIR indexation step, using the input fof:

```bash
./StrainFLAIR.sh index -i data/list_fasta_sim.txt -o myproject
```

The reference graph and its associated additionnal files have been created in the new directory `myproject`.

**Third**, query the variation graph by mapping any reads on it. Reads can be compressed in a fastq.gz format. 

```bash
./StrainFLAIR.sh query -g myproject/graphs/all_graphs -f1 data/manuscript_simulated_dataset/mix_o104_iai39_k12_100k_e0.001.fastq.gz -t 24 -p myproject/graphs/dict_clusters.pickle -d myproject -o manuscriptsim
```

That's it! Mapping files are available in `mapping/`, and abundance tables are available in `results/`.

The final result is a csv table containing each reference genome in line identified by their accession number. Columns contained the proportion of detected genes and the estimated strain-level abundance according to different computation methods. Here we found close estimations compared to what is expected: 44.7% for O104:H4 2011C-3493 (NC_018658.1), 35% for IAI39 (NC_011750.1), 20.3% for K-12 MG1655, and 0 for all the absent strains.

Here is the complete output:
```
,detected_genes,mean_abund,mean_abund_nz,median_abund,median_abund_nz
NC_000913.3,0.9787234042553191,20.29226977437424,21.05505540046613,19.25596456585335,19.939248305563503
NC_002695.2,0.2019113460756405,0.0,0.0,0.0,0.0
NC_018658.1,0.9935683046050939,44.67232236815524,44.1192172117696,46.3219304819112,45.76611823931458
NC_011750.1,0.9885301614273577,35.03540785747052,34.825727387764275,34.42210495223545,34.29463345512191
NZ_CP028116.1,0.22883087400681043,0.0,0.0,0.0,0.0
NZ_CP007592.1,0.20748792270531402,0.0,0.0,0.0,0.0
NC_013654.1,0.15242438790206433,0.0,0.0,0.0,0.0
```


**Full commands:**

```
find data/manuscript_simulated_dataset/references/* -not -name "*bl21*" > data/list_fasta_sim.txt
./StrainFLAIR.sh index -i data/list_fasta_sim.txt -o myproject
./StrainFLAIR.sh query -g myproject/graphs/all_graphs -f1 data/manuscript_simulated_dataset/mix_o104_iai39_k12_100k_e0.001.fastq.gz -t 24 -p myproject/graphs/dict_clusters.pickle -d myproject -o manuscriptsim
```
