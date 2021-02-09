# StrainFLAIR cookbook

- - - -
**Prerequisite**

Datasets for this cookbook are available in the folder `data/`. It contains 4 complete genomes (fasta files) and a mixture of simulated reads with 0.001% of errors composed of 67% of reads from D4 and 33% of reads from LM33 (fastq file).

**Full commands:** For each example, the full commands are proposed at the end of the section.

**Computation time** is approximatel 3 minutes for the indexation step, 3 minutes for the mapping step and 1 minute for the query step.
- - - -

Consider 4 reference genomes:

**First**, create a *file of file* (fof) used as input for StrainFLAIR. Each line of the input fof contains exactly one fasta name. We may create this fof as follows:

```bash
ls data/*.fasta > list_fasta.txt
```

**Second**, run StrainFLAIR indexation step, using the input fof:

```bash
./StrainFLAIR.sh index -i list_fasta.txt -o myproject
```

The reference graph and its associated additionnal files have been created in the new directory `myproject`.

**Third**, map any reads on the so created reference graph. Reads can be compressed in a fastq.gz format. 
The mapping output needs to be converted into a JSON file using `vg view`.

```bash
vg mpmap -x myproject/graphs/all_graphs.xg -g myproject/graphs/all_graphs.gcsa -s myproject/graphs/all_graphs.snarls -f myreads.fastq -t 24 -M 10 -m -L 0 > myproject/mapping_output.gamp
vg view -j -K myproject/mapping_output.gamp > myproject/mapping_output.json
```

**Fourth**, query the mapped reads to the pangenome graph.

```bash
./StrainFLAIR.sh query -g myproject/graphs/all_graphs.gfa -m myproject/mapping_output.json -p test/graphs/dict_clusters.pickle -o myproject/output
```

That's it! Results are available in `output/`.

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
./StrainFLAIR.sh index -i list_fasta.txt -o myproject

vg mpmap -x myproject/graphs/all_graphs.xg -g myproject/graphs/all_graphs.gcsa -s myproject/graphs/all_graphs.snarls -f myreads.fastq -t 24 -M 10 -m -L 0 > myproject/mapping_output.gamp
vg view -j -K myproject/mapping_output.gamp > myproject/mapping_output.json

./StrainFLAIR.sh query -g myproject/graphs/all_graphs.gfa -m myproject/mapping_output.json -p test/graphs/dict_clusters.pickle -o myproject/output
```
