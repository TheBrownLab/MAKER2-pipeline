########################################################################################################################
# Robert Jones (rej110@msstate.edu)
# usage: maker_run.py {pass #} {# of threads}
#
# example: maker_run.py 2 10
########################################################################################################################
import argparse
import os
import re
import shutil
import subprocess


class Paths:
    def __init__(self, pass_num):
        # i = iteration; p = Pass Number
        self.base = '/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Maker_alt_EST/ref10'

        # Info about genome here
        self.base_name = f'oneContg.ref10'
        self.genome_path = f'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/PseudoContig/oneContg.ref10.fa'

        # RNA Seq
        self.rna_seq = f'\
        /mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Acan_NEFF/data/sra_assembly/RNAseq_trinity_dal/Trinity.fasta'

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


def bash_command(arg):
    subprocess.call(f'{arg}', shell=True, executable='/bin/bash')


# This fucntion generates the RepeatModeler BLAST databases to be used in RepeatModeler
def gen_repmod_blastdb(exes, my_paths, i):
    os.chdir(my_paths.rep_mod_db)

    shutil.copyfile(my_paths.genome_path, f'./db{i}.fa')

    # Use absolute path to BuildDatabase executable. Change if not run on lugh.
    exes.build_db(arg=f'-name db{i} {my_paths.genome_path}')



def make_dir(dir_path):
    try:
        os.mkdir(dir_path)
    except FileExistsError:
        shutil.rmtree(dir_path)
        os.mkdir(dir_path)


def get_aug_species(i, pass_num):
    my_paths = Paths(pass_num=pass_num)

    busco_path = my_paths.busco_dir
    aug_config = my_paths.aug_config_path
    species = my_paths.base_name

    make_dir(f'{aug_config}/species/{species}')
    print(f'{aug_config}/species/{species}')

    endings = ['_intron_probs.pbl', '_exon_probs.pbl', '_parameters.cfg.orig1', '_igenic_probs.pbl',
               '_parameters.cfg',
               '_weightmatrix.txt', '_metapars.cfg', '_metapars.utr.cfg', '_metapars.cgp.cfg']

    # r=root, d=directories, f = files
    for root, dirs, files in os.walk(busco_path):

        for file in files:
            for ending in endings:
                new_file = f'{aug_config}/species/{species}/{species}{ending}'

                if file.endswith(ending) is True:
                    file = os.path.join(root, file)
                    shutil.copy(file, new_file)
                    if ending is '_parameters.cfg.orig1' or ending is '_parameters.cfg':
                        old_prefix = file.split('/')[-1]
                        old_prefix = old_prefix.split('_parameters')[0]

                        f = open(new_file, 'r')
                        filedata = f.read()
                        f.close()

                        newdata = filedata.replace(old_prefix, species)

                        f = open(new_file, 'w')
                        f.write(newdata)
                        f.close()

    return species


def run_genemark(threads):
    my_paths = Paths(pass_num=1)
    make_dir(my_paths.gm_dir)
    os.chdir(my_paths.gm_dir)

    cmd = f'''\
        conda activate genemark_es
        
        gmes_petap.pl --ES --min_contig 10000 --cores {threads} --sequence {my_paths.genome_path} \
        '''
    bash_command(cmd)


def run_snap():
    # Commands to train SNAP
    # these lines go in the job script that follows if pass 1 or 2
    make_dir(new_paths.snap_dir)
    cmd = f'''        
            cd {new_paths.snap_dir}
            maker2zff {new_paths.maker_all_gff}
            fathom -categorize 1000 genome.ann genome.dna
            fathom -export 1000 -plus uni.ann uni.dna 
            forge export.ann export.dna
            hmm-assembler.pl {new_paths.base_name}_snap . > {new_paths.base_name}_snap.hmm
        '''
    bash_command(cmd)


def run_tran_busco(threads):
    # BUSCO run after passes 1 and 2
    # Trains AUGUSTUS
    make_dir(new_paths.busco_dir)
    cmd = f'''
            cd {new_paths.busco_dir}
            export AUGUSTUS_CONFIG_PATH={new_paths.aug_config_path}
            python {my_exes.busco_exe} \\
                -i {new_paths.maker_tran} \\
                -l {new_paths.busco_lineage} \\
                -o {new_paths.base_name}_tran \\
                -m genome \\
                -c {threads} \\
                --long \\
                --augustus_parameters='--progress=true'
        '''
    bash_command(cmd)


def run_prot_busco(threads):
    # Final BUSCO run (likely needs to change to pass_num == 3)
    # fourth pass through maker likely not necesarry
    # Does NOT train AUGUSTUS
    cmd = f'''
            cd {new_paths.final_busco_dir}
            export AUGUSTUS_CONFIG_PATH={new_paths.aug_config_path}
            python {my_exes.busco_exe} \\
                -i {new_paths.maker_prot} \\
                -l {new_paths.busco_lineage} \\
                -o {new_paths.base_name}_prot \\
                -m proteins \\
                -c {threads}
        '''
    bash_command(cmd)


def run_maker(paths, threads):
    # Makes MAKER directory and cd's to it
    make_dir(new_paths.maker_dir)
    os.chdir(new_paths.maker_dir)

    # Generates and edits MAKER ctl files
    maker_ctl()

    # Runs Maker
    cmd = f'''\
        cd {paths.maker_dir}
        maker -c {threads - 1} -fix_nucleotides maker_opts_edit.ctl maker_bopts.ctl mod_maker_exe.ctl

        cd {paths.maker_out}
        fasta_merge -d {paths.maker_log}
        gff3_merge -d {paths.maker_log}                            
        '''
    bash_command(cmd)


def mod_maker_exe():
    with open('mod_maker_exe.ctl', 'w') as new_exe_file, open('maker_exe.ctl', 'r') as exe_file:
        for line in exe_file:
            line = line.strip()
            if line.split('=')[0] == 'RepeatMasker':
                line = line.split('=')[0] + '=' + my_exes.rep_mask_exe + ' #' + line.split('#')[1]
            new_exe_file.write(line + '\n')


def maker_ctl(pass_num, est_alt):
    # Create the 3 MAKER2 .ctl files
    my_exes.maker(arg='-CTL')
    mod_maker_exe()

    # Edit maker_opts.ctl - This file contatins the options for the maker run.
    with open(f'maker_opts_edit.ctl', 'w') as new_opt_file, open('maker_opts.ctl', 'r') as opt_file:
        if est_alt is False:
            est_key = 'est'
        elif est_alt is True:
            est_key = 'altest'

        all_passes = {'genome': old_paths.genome_path, 'est_pass': 1, 'rmlib': old_paths.rep_mod_out,
                      'gmhmm': old_paths.gm_hmm, est_key: old_paths.rna_seq}
        pass_1_key = {'est2genome': 1, 'protein2genome': 1, 'keep_preds': 1, 'single_exon': 1}
        pass_2_key = {'maker_gff': old_paths.maker_all_gff, 'snaphmm': old_paths.snap_hmm,
                      'augustus_species': old_paths.base_name}

        for line in opt_file:
            line = line.strip()
            line_begin = line.split('=')[0]
            try:
                line_end = line.split('#')[1]
            except IndexError:
                pass

            if line_begin in all_passes.keys():
                line = f'{line_begin}={all_passes[line_begin]} #{line_end}'
            elif line_begin in pass_1_key.keys() and pass_num == 1:
                line = f'{line_begin}={pass_1_key[line_begin]} #{line_end}'
            elif line_begin in pass_1_key.keys() and pass_num == 2:
                line = f'{line_begin}={pass_2_key[line_begin]} #{line_end}'

            new_opt_file.write(line + '\n')

    # Removes orginal, unedited maker_opts.ctl and maker_exe.ctl
    os.remove('maker_opts.ctl')
    os.remove('maker_exe.ctl')


def sumamry():

    dirs = ['BUSCO_Pass1/PseudoContigs/oneContg.', 'BUSCO_Pass2/PseudoContigs/oneContg.', 'Final_BUSCO/']
    base = '/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/MakerPipeline'
    error1 = []
    error2 = []
    error3 = []
    with open('BUSCO_Scores.csv', 'w') as outfile:
        for i in range(1, 23):
            complete_list = []
            for x in dirs:
                if x is not 'Final_BUSCO/':
                    file = f'{base}/{x}ref{i}/run_{Sets.key[i]}/short_summary_{Sets.key[i]}.txt'
                else:
                    file = f'{base}/{x}ref{i}/run_oneContg.ref{i}_BUSCO/short_summary_oneContg.ref{i}_BUSCO.txt'

                with open(file, 'r') as infile:
                    for line in infile:
                        line = line.strip()

                        complete = re.findall(r'C:[0-9]{0,2}\.[0-9]', line)
                        if complete:
                            complete_list.append(complete[0].split(':')[1])

            if float(complete_list[0]) == 0:
                error1.append(int(i))
            if float(complete_list[1]) == 0:
                error2.append(int(i))
            if float(complete_list[2]) == 0:
                error3.append(int(i))

            outfile.write(f'ref{i},{key[i]},{complete_list[0]},{complete_list[1]},{complete_list[2]}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for running MAKER2 Pipeline.',
                                     usage="maker_run.py -p PASSAGE -t THREADS")
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads, default:1')
    parser.add_argument('-p', '--passage', type=int, default=1,
                        help='Passage number through MAKER2 pipeline (i.e. 1, 2, or 3)')
    parser.add_argument('-e', '--est_alt', type=bool, default=False,
                        help='Passage number through MAKER2 pipeline (i.e. 1, 2, or 3)')
    args = parser.parse_args()

    new_paths = Paths(pass_num=args.passage)
    old_paths = Paths(pass_num=(args.passage - 1))
    my_exes = Executables()

    get_aug_species(pass_num=args.passage - 1)
    run_maker(paths=new_paths, threads=args.threads)

    if args.passage < 3:
        run_snap()

        run_tran_busco(threads=args.threads)

    if args.passage == 3:
        run_prot_busco(threads=args.threads)


