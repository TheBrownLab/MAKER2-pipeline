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
`wget https://github.com/rej110/myMAKER2-Pipeline.git` 

## 2. Usage
##### 2.1 Data
* Place genome, est data, and protein homology data in `data/` directory
    * Name genome as `**orgname**.genome.fas`
    * Name est data as `**orgname**.est.fas`
    * Name protein homology data as `**orgname**.protein.fas`
* Where `**orgname**` is the name of the organism. This will be different if est and protein data originates from alternate organism
* See `example/data/` to see where and how data files should be placed and named.

##### 2.2 Command Line
* `python3 maker_run.py [options]`
    * `-p PASSAGE, --passage PASSAGE` : Passage number through MAKER2 pipeline (i.e. 1, 2, or 3), default:1
    * `-t THREADS, --threads THREADS` : Number of threads, default:1
    * `-a, --alt_est` : If est and protein homology data come from an alternate organims, default:False

##### 2.3 Work Flow
* Pass 1 - Evidence based predicitons only
    * Receives:
        * EST data
        * Protein homology data

        
* Pass 2 - _Ab initio_ predictions only.
    * Receives:
        * GeneMarkES Training
        * SNAP Training
        * Augustus Species 

 
* Pass 3 - _Ab initio_ predictions only.
    * Receives:
        * GeneMarkES Training
        * SNAP Training
        * Augustus Species 
    
    
* Pass 4 - Evidence based and _Ab initio_ predictions
    * Receives:
        * GeneMarkES Training
        * SNAP Training
        * Augustus Species
        * EST data
        * Protein homology data


SNAP and Augustus are retrained between passes\
GeneMarkES is self trained
### References
