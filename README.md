# myMAKER2-Pipeline

## 1. Installation
##### 1.1 Install Conda
* Install Conda via Miniconda or Anaconda: \
https://docs.conda.io/projects/conda/en/latest/user-guide/install/
* Prepare conda virtual enviornment: \
`conda env create -f maker_env.yml`
* Activate maker environment: \
`conda activate maker`
##### 1.2  Install myMAKER2-Pipeline
* Download from GitHub: \
`wget https://github.com/rej110/myMAKER2-Pipeline.git` 

## 2. Usage
##### 2.1 Data
* Place genome, est data, and protein homology data in `data/` directory
    * Name genome as `**orgname**.genome.fas`
    * Name est data as `**orgname**.est.fas`
    * Name protein homology data as `**orgname**.protein.fas`
        * Where `**orgname**` is the name of the organism
        * See `example/data/` to see where and how data files should be placed and named.

##### 2.2 Command Line
* `python3 maker_run.py [options]`
    * `-p PASSAGE, --passage PASSAGE` : Passage number through MAKER2 pipeline (i.e. 1, 2, or 3), default:1
    * `-t THREADS, --threads THREADS` : Number of threads, default:1
    * `-e, --est_alt` : est and protein homology data dervived from an alternate organims


### References
