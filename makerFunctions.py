########################################################################################################################
#
#
#
########################################################################################################################

import re
import os
import sys
import shutil
import subprocess

# sys.path.append(r'/mnt/home/rej110/Scripts/MakerPipeline')
import makerSettings as Sets


def bash_command(arg):
    subprocess.call(f'{arg}', shell=True, executable='/bin/bash')


# This fucntion generates the RepeatModeler BLAST databases to be used in RepeatModeler
def gen_repmod_blastdb(my_exes, my_paths, i):
    os.chdir(my_paths.rep_mod_db)

    shutil.copyfile(my_paths.genome_path, f'./db{i}.fa')

    # Use absolute path to BuildDatabase executable. Change if not run on lugh.
    my_exes.build_db(arg=f'-name db{i} {my_paths.genome_path}')


def make_dir(dir_path):
    try:
        os.mkdir(dir_path)
    except FileExistsError:
        shutil.rmtree(dir_path)
        os.mkdir(dir_path)


def get_aug_species(i, pass_num):
    my_paths = Sets.Paths(i=i, pass_num=pass_num)

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
