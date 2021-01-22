## BUILD PACKAGE FROM SOURCE:
```bash
git clone https://github.com/kevsilva/StrainFLAIR
conda create -p strainflair-env
conda activate strainflair-env
conda install conda-build -c bioconda 
conda-build ./StrainFLAIR/conda/strainflair -c bioconda
```
