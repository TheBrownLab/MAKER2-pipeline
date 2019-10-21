########################################################################################################################
#
#
#
########################################################################################################################
import os
import sys
#               Path to MakerPipline modules (makerFunctions and makerSettings)
# sys.path.append(r'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Acan_NEFF/NEFF_maker/MakerScripts')
import makerFunctions as Funcs
import makerSettings as Sets


def main():
    threads = int(sys.argv[1])

    for i in [1]:
        my_paths = Sets.Paths(i=i, pass_num=1)
        my_exes = Sets.Executables()

        # Creates GeneMark directory, build BLAST database, cd to dir that was just made
        Funcs.make_dir(my_paths.rep_mod_db)
        Funcs.gen_repmod_blastdb(my_paths=my_paths, my_exes=my_exes, i=i)
        os.chdir(my_paths.rep_mod_db)

        with open(f'RMod{i}.sh', 'w') as job_script:
            job_script.write(f'''\
#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threads {threads}

module load RepeatModeler

cd {my_paths.rep_mod_db}                
RepeatModeler -pa {threads-1} -database db{i}
\
''')

        # Submit job script
        my_exes.qsub(arg=f'RMod{i}.sh')
########################################################################################################################


main()
