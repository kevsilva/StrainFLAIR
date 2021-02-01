# StrainFLAIR cookbook

TODO

- - - -
**Prerequisite**

Datasets for this cookbook are available in the folder TODO

**Full commands:** For each example, the full commands are proposed at the end of the section.

**Computation time** TODO.
- - - -

Consider 4 reference genomes :

**First**, creating a *file of file* (fof) is needed, used as input for StrainFLAIR. Each line of the input fof contains exactly one fasta name. We may create this fof as following:

```bash
ls *.fasta > list_fasta.txt
```

**Second**, run StrainFLAIR indexation step, using the input fof:

```bash
./StrainFLAIR.sh index -i list_fasta.txt -o myproject
```

The pangenome graph and its required additionnal files have been created.

**Third**, map your reads onto the pangenome graph (can be compressed in a fastq.gz format). The mapping output needs to be converted into a JSON file.

```bash
vg mpmap -x myproject/graphs/all_graphs.xg -g myproject/graphs/all_graphs.gcsa -s myproject/graphs/all_graphs.snarls -f myreads.fastq -t 24 -M 10 -m -L 0 > myproject/mapping_output.gamp
vg view -j -K myproject/mapping_output.gamp > myproject/mapping_output.json
```

**Fourth**, query the mapped reads to the pangenome graph.

```bash
./StrainFLAIR.sh query -g myproject/graphs/all_graphs.gfa -m myproject/mapping_output.json -p test/graphs/dict_clusters.pickle -o myproject/output
```

That's it! Results are available in `output/`.

(ouput example)

**Full commands:**

```
./StrainFLAIR.sh index -i list_fasta.txt -o myproject

vg mpmap -x myproject/graphs/all_graphs.xg -g myproject/graphs/all_graphs.gcsa -s myproject/graphs/all_graphs.snarls -f myreads.fastq -t 24 -M 10 -m -L 0 > myproject/mapping_output.gamp
vg view -j -K myproject/mapping_output.gamp > myproject/mapping_output.json

./StrainFLAIR.sh query -g myproject/graphs/all_graphs.gfa -m myproject/mapping_output.json -p test/graphs/dict_clusters.pickle -o myproject/output
```
