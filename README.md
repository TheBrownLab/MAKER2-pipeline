# myMAKER2-Pipeline

## 1. Installation
##### 1.1 Install Dependencies 
* RepeatModeler
* RepeatMasker
* MAKER2
* GenemarkES
* Augustus
* SNAP
* BUSCO
##### 1.2  Install myMAKER2-Pipeline
* Download from GitHub: \


## 2. Usage
#### 2.1 Data
* Place genome, est data, and protein homology data in `data/` directory
    * Name genome as `**orgname**.genome.fas`
    * Name est data as `**orgname**.est.fas`
    * Name protein homology data as `**orgname**.protein.fas`
* Where `**orgname**` is the name of the organism. This will be different if est and protein data originates from alternate organism
* See `example/data/` to see where and how data files should be placed and named.

#### 2.2 Command Line
* `python3 maker_run.py [options]`
    * `-p PASSAGE, --passage PASSAGE` : Passage number through MAKER2 pipeline (i.e. 1, 2, or 3), default:1
    * `-t THREADS, --threads THREADS` : Number of threads, default:1
    * `-a, --alt_est` : If est and protein homology data come from an alternate organims, default:False

#### 2.3 MAKER2 ctl
* All Passes
    * genome = _path to genome_
    * rm_lib = _path to RepeatModeler output_
    * protein = _path to protein homology data_
    * est = path to est data from same organism (blank if from alternate organism)
    * altest = path to est data from alternate organism (blank if from same organism)
    
    
* Passes 2, 3, & 4
    * snaphmm = path to SNAP hmm
    * gmhmm = path to GeneMark hmm
    * augustus_species = name of AUGUSTUS species (species model produced by BUSCO)
    * maker_gff = path to MAKER2 derived gff
    
    
* Pass 1 (only evidence based gene models are reported)
    * est2genome = 1
    * protein2genome = 1
    * keep_preps = 0
    
    
* Passes 2 & 3 (only gene models supported by homology evidence are reported)
    * est2genome = 0
    * protein2genome = 0
    * keep_preps = 0
    
    
* Pass 4 (all gene models are reported)
    * est2genome = 0
    * protein2genome = 0
    * keep_preps = 1

*SNAP and Augustus are retrained between passes\
*GeneMarkES is self trained
### References
