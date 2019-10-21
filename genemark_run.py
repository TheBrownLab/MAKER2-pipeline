########################################################################################################################
#
#
#
########################################################################################################################
import os
import sys
import inspect

# sys.path.append(r'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Acan_NEFF/NEFF_maker/MakerScripts')
import makerFunctions as Funcs
import makerSettings as Sets


def main():
    threads = int(sys.argv[1])

    for i in [1]:
        my_paths = Sets.Paths(i=i, pass_num=1)
        my_exes = Sets.Executables()

        Funcs.make_dir(my_paths.gm_dir)

        os.chdir(my_paths.gm_dir)

        with open(f'GeneMark{i}.sh', 'w') as job_script:
            job_str = inspect.cleandoc(f'''\
                #!/bin/bash
                #$ -S /bin/bash
                . /etc/profile
                #$ -cwd
                #$ -pe threads {threads}

                conda activate genemark_es

                gmes_petap.pl --ES --min_contig 10000 --cores 20 --sequence {my_paths.genome_path} \
                            ''')
            job_script.write(job_str)

        # Submit job script
        my_exes.qsub(f'GeneMark{i}.sh')


main()
