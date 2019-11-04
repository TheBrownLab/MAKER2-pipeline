########################################################################################################################
# Robert Jones (rej110@msstate.edu)
# usage: maker_run.py [options]
#        run maker_run.py -h for options
########################################################################################################################
import os
import re
import sys
import shutil
import argparse
import subprocess


class Paths:
    # Executables
    rep_mask_exe = '/mnt/home/software/sources/RepeatMasker/4-0-8/RepeatMasker'

    def __init__(self, pass_num):
        self.base = os.getcwd()
        # Gets BaseName and paths to genome, est, and protein homology data
        self.base_name, self.genome, self.est, self.prot = self.get_data_paths()

        # RepeatModeler
        self.rep_mod_db = f'{self.base}/RepeatModeler'
        self.rep_mod_out = f'{self.rep_mod_db}/{self.base_name}.db-families.fa'

        # GeneMarkES
        self.gm_dir = f'{self.base}/GeneMarkES/'
        self.gm_hmm = f'{self.gm_dir}/output/gmhmm.mod'

        # Path of current Maker run
        self.maker_dir = f'{self.base}/Maker_Pass{pass_num}'
        self.maker_out = f'{self.maker_dir}/{self.base_name}.maker.output'
        self.maker_log = f'{self.maker_out}/{self.base_name}_master_datastore_index.log'
        self.maker_all_gff = f'{self.maker_out}/{self.base_name}.all.gff'
        self.maker_tran = f'{self.maker_out}/{self.base_name}.all.maker.transcripts.fasta'
        self.maker_prot = f'{self.maker_out}/{self.base_name}.all.maker.proteins.fasta'

        # Path to SNAP training
        self.snap_dir = f'{self.base}/Snap_Training{pass_num}'
        self.snap_hmm = f'{self.snap_dir}/{self.base_name}_snap.hmm'

        # Path to files for BUSCO/AugustusTraining
        self.busco_lineage = f'/mnt/scratch/brownlab/rej110/.conda/envs/maker/eukaryota_odb9'
        self.busco_dir = f'{self.base}/BUSCO_Pass{pass_num}/'
        self.aug_config_path = '/mnt/scratch/brownlab/rej110/.conda/envs/maker/config'

    def get_data_paths(self):
        base, genome, est, prot = None, None, None, None
        for root, dirs, files in os.walk(self.base):
            for filename in files:
                if filename.endswith('.genome.fas'):
                    base = os.path.splitext(filename)[0]  # .split('.')[0]
                    genome = os.path.join(root, filename)
                elif filename.endswith('.est.fas'):
                    est = os.path.join(root, filename)
                elif filename.endswith('protein.fas'):
                    prot = os.path.join(root, filename)

        if base is None or genome is None or est is None or prot is None:
            sys.exit('Data file(s) are missing')
        else:
            return base, genome, est, prot


def execute(cmd):
    subprocess.run(cmd, shell=True, executable='/bin/bash')


def run_rep_mod(threads, new):
    make_dir(new.rep_mod_db)
    os.chdir(new.rep_mod_db)

    shutil.copyfile(new.genome, f'./{new.base_name}.db.fas')

    cmd = f'''
    module load RepeatModeler/1.0.11
    BuildDatabase -name {new.base_name}.db ./{new.base_name}.db.fas
    RepeatModeler -pa {threads - 1} -database {new.base_name}.db
    module unload RepeatModeler/1.0.11
    '''
    execute(cmd)


def make_dir(dir_path):
    try:
        os.mkdir(dir_path)
    except FileExistsError:
        shutil.rmtree(dir_path)
        os.mkdir(dir_path)


def get_aug_species(new):

    busco_path = new.busco_dir
    aug_config = new.aug_config_path
    species = new.base_name

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

    cmd = f'''
    gmes_petap.pl --ES --min_contig 10000 --cores {threads} --sequence {my_paths.genome}
    '''
    execute(cmd)


def run_snap(new):
    # Commands to train SNAP
    # these lines go in the job script that follows if pass 1 or 2
    make_dir(new.snap_dir)
    os.chdir(new.snap_dir)

    cmd = f'''
    maker2zff {new.maker_all_gff}
    fathom -categorize 1000 genome.ann genome.dna
    fathom -export 1000 -plus uni.ann uni.dna
    forge export.ann export.dna
    hmm-assembler.pl {new.base_name}_snap . > {new.base_name}_snap.hmm
    '''
    execute(cmd)


def run_busco(threads, new):
    # Final BUSCO run (likely needs to change to pass_num == 3)
    # fourth pass through maker likely not necesarry
    # Does NOT train AUGUSTUS
    make_dir(new.busco_dir)
    os.chdir(new.busco_dir)

    cmd = f'''
    export AUGUSTUS_CONFIG_PATH={new.aug_config_path}
    run_BUSCO.py \\
    -i {new.maker_prot}  \\
    -l {new.busco_lineage} \\
    -o {new.base_name}_prot \\
    -m proteins \\
    -c {threads}
    
    run_BUSCO.py \\
    -i {new.maker_tran} \\
    -l {new.busco_lineage} \\
    -o {new.base_name}_tran  \\
    -m genome \\
     -c {threads} \\
    --long
    '''
    execute(cmd)


def run_maker(new, old, threads, pass_num, alt_est):
    # Makes MAKER directory and cd's to it
    make_dir(new.maker_dir)
    os.chdir(new.maker_dir)

    # Generates and edits MAKER ctl files
    maker_ctl(pass_num=pass_num, alt_est=alt_est, old=old)

    # Runs Maker
    cmd = f'''
    maker -c {threads - 1} -fix_nucleotides maker_opts_edit.ctl maker_bopts.ctl mod_maker_exe.ctl
    cd {new.maker_out}
    fasta_merge -d {new.maker_log}
    gff3_merge -d {new.maker_log}
    '''

    execute(cmd)


def mod_maker_exe(old):
    with open('mod_maker_exe.ctl', 'w') as new_exe_file, open('maker_exe.ctl', 'r') as exe_file:
        for line in exe_file:
            line = line.strip()
            if line.split('=')[0] == 'RepeatMasker':
                line = line.split('=')[0] + '=' + old.rep_mask_exe + ' #' + line.split('#')[1]
            new_exe_file.write(line + '\n')
    os.remove('maker_exe.ctl')


def maker_ctl(pass_num, alt_est, old):
    # Create the 3 MAKER2 .ctl fils
    cmd = 'maker -CTL'
    execute(cmd)

    mod_maker_exe(old=old)

    # Sets paths to genome, est evidence, and protein homology evidence in maker ctl file.
    # For the first pass of maker. Creates gene model needed to train SNAP and AUGUSTUS
    all_passes = {'genome': old.genome, 'rmlib': old.rep_mod_out, 'protein': old.prot}
    # Sets est evidence to est or altest depending if est evidence is from an alternate organism
    if alt_est is False:
        all_passes['est'] = old.est
    elif alt_est is True:
        all_passes['altest'] = old.est

    # Sets paths to SNAP hmm, Genemark hmm, and the AUGUSTUS species model.
    # For the "true passes" of maker. Passes 2, 3, and 4
    pass_234_key = {'snaphmm': old.snap_hmm, 'gmhmm': old.gm_hmm, 'augustus_species': old.base_name,
                    'maker_gff': old.maker_all_gff}

    pass_1_key = {'est2genome': 1, 'protein2genome': 1, 'keep_preds': 0}
    pass_23_key = {'est2genome': 0, 'protein2genome': 0, 'keep_preds': 0}
    pass_4_key = {'est2genome': 0, 'protein2genome': 0, 'keep_preds': 1}

    # Edit maker_opts.ctl - This file contatins the options for the maker run
    with open(f'maker_opts_edit.ctl', 'w') as new_opt_file, open('maker_opts.ctl', 'r') as opt_file:
        for line in opt_file:
            line = line.strip()
            line_begin = line.split('=')[0]
            try:
                line_end = line.split('#')[1]
            except IndexError:
                pass

            try:
                if line_begin in all_passes.keys():
                    line = f'{line_begin}={all_passes[line_begin]} #{line_end}'

                elif line_begin in pass_1_key.keys() and pass_num == 1:
                    line = f'{line_begin}={pass_1_key[line_begin]} #{line_end}'

                elif line_begin in pass_234_key.keys() and (pass_num == 2 or pass_num == 3 or pass_num == 4):
                    line = f'{line_begin}={pass_234_key[line_begin]} #{line_end}'

                elif line_begin in pass_23_key.keys() and (pass_num == 2 or pass_num == 3):
                    line = f'{line_begin}={pass_23_key[line_begin]} #{line_end}'

                elif line_begin in pass_4_key.keys() and pass_num == 4:
                    line = f'{line_begin}={pass_4_key[line_begin]} #{line_end}'

            except KeyError:
                pass

            new_opt_file.write(line + '\n')

    # Removes orginal, unedited maker_opts.ctl and maker_exe.ctl
    os.remove('maker_opts.ctl')


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
                        help='Passage through MAKER2 pipeline to start from (i.e. 1, 2, 3, or 4; default=1)')
    parser.add_argument('-a', '--alt_est', default=False, action='store_true',
                        help='If est and protein data from an alternate organism')
    args = parser.parse_args()
    passage = args.passage

    for this_pass in range(passage, 5):
        new_paths = Paths(pass_num=this_pass)
        old_paths = Paths(pass_num=(this_pass - 1))

        if this_pass == 1:
            # run_genemark(threads=args.threads)
            run_rep_mod(threads=args.threads, new=new_paths)
            run_maker(threads=args.threads, pass_num=this_pass, alt_est=args.alt_est,
                      new=new_paths, old=old_paths)
            run_snap(new=new_paths)
            run_busco(threads=args.threads, new=new_paths)
            get_aug_species(new=new_paths)
            os.chdir(new_paths.base)

        elif this_pass == 2 or this_pass == 3:
            run_maker(threads=args.threads, pass_num=this_pass, alt_est=args.alt_est,
                      new=new_paths, old=old_paths)
            run_snap(new=new_paths)
            run_busco(threads=args.threads, new=new_paths)
            get_aug_species(new=new_paths)
            os.chdir(new_paths.base)

        elif this_pass == 4:
            run_maker(threads=args.threads, pass_num=this_pass, alt_est=args.alt_est,
                      new=new_paths, old=old_paths)
            run_busco(threads=args.threads, new=new_paths)
            os.chdir(new_paths.base)
