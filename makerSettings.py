########################################################################################################################
#
#
#
######################################################################################################################
import subprocess
import os


########################################################################################################################
class Executables:
    def __init__(self):
        # Paths to executables executed outside of these python scripts (i.e. jobs scripts, control files)
        self.rep_mask_exe = '/mnt/home/software/sources/RepeatMasker/4-0-8/RepeatMasker'
        self.busco_exe = '/mnt/scratch/brownlab/rej110/.conda/envs/maker/bin/run_BUSCO.py'
        self.bedtool_exe = '/mnt/scratch/brownlab/rej110/.conda/envs/maker/bin/bedtools'

        # Paths to executables executed within these python scripts
        # Functions to execute follow
        self.build_db_exe = '/mnt/home/software/sources/RepeatModeler-open-1.0.11//BuildDatabase'
        self.maker_exe = '/mnt/home/software/sources/maker/bin/maker'
        self.qsub_exe = '/opt/sge/bin/lx-amd64//qsub'
        self.fasta_merge_exe = '/mnt/home/software/anaconda/anaconda2/bin/fasta_merge'
        self.gff3_merge_exe = '/mnt/home/software/anaconda/anaconda2/bin/gff3_merge'
        self.maker2zff_exe = '/mnt/home/software/anaconda/anaconda2/bin/maker2zff'
        self.make_aug_model = '/mnt/home/rej110/Scripts/MakerPipeline/MakerPipeline/make_AugustusModel.py'

    # Functions to execute above executables
    # Do Not Modify
    def build_db(self, arg):
        subprocess.call(f'{self.build_db_exe} {arg}', shell=True, executable='/bin/bash')

    def qsub(self, arg):
        subprocess.call(f'{self.qsub_exe} {arg}', shell=True, executable='/bin/bash')

    def maker(self, arg):
        subprocess.call(f'{self.maker_exe} {arg}', shell=True, executable='/bin/bash')

    def fasta_merge(self, arg):
        subprocess.call(f'{self.fasta_merge_exe} {arg}', shell=True, executable='/bin/bash')

    def gff3_merge(self, arg):
        subprocess.call(f'{self.gff3_merge_exe} {arg}', shell=True, executable='/bin/bash')

    def maker2zff(self, arg):
        subprocess.call(f'{self.maker2zff_exe} {arg}', shell=True, executable='/bin/bash')


########################################################################################################################


########################################################################################################################
class Paths:
    def __init__(self, i, pass_num):
        # i = iteration; p = Pass Number

        self.base = '/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Maker_alt_EST/ref10'

        # Info about genome here
        self.base_name = f'oneContg.ref10'
        self.genome_path = f'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/PseudoContig/oneContg.ref10.fa'

        # RNA Seq
        self.rna_seq = f'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Acan_NEFF/data/sra_assembly/RNAseq_trinity_dal/Trinity.fasta'

        # RepeatModeler
        self.rep_mod_db = f'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/MakerPipeline/RepeatModeler/ref10'
        self.rep_mod_out = f'{self.rep_mod_db}/db10-families.fa'

        # GeneMarkES
        self.gm_dir = f'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/MakerPipeline/GeneMarkES/ref10'
        self.gm_hmm = f'{self.gm_dir}/output/gmhmm.mod'

        # Path of current Maker run
        self.maker_dir = f'{self.base}/maker/Maker_Pass{pass_num}'
        self.maker_out = f'{self.maker_dir}/{self.base_name}.maker.output'
        self.maker_log = f'{self.maker_out}/{self.base_name}_master_datastore_index.log'
        self.maker_all_gff = f'{self.maker_out}/{self.base_name}.all.gff'
        self.maker_tran = f'{self.maker_out}/{self.base_name}.all.maker.transcripts.fasta'
        self.maker_prot = f'{self.maker_out}/{self.base_name}.all.maker.proteins.fasta'

        # Path to SNAP training
        self.snap_dir = f'{self.base}/maker/Snap_Training{pass_num}'
        self.snap_hmm = f'{self.snap_dir}/ref10_snap.hmm'

        # Path to files for BUSCO/AugustusTraining
        self.busco_lineage = f'/mnt/scratch/brownlab/rej110/.conda/envs/maker/eukaryota_odb9'
        self.busco_dir = f'{self.base}/maker/BUSCO_Pass{pass_num}/'
        self.aug_config_path = '/mnt/scratch/brownlab/rej110/.conda/envs/maker/config'

        # Path to Post Maker files
        self.final_busco_dir = f'{self.base}/maker/Final_BUSCO/'


########################################################################################################################
