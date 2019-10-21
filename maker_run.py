########################################################################################################################
# Robert Jones (rej110@msstate.edu)
# usage: maker_run.py {pass #} {# of threads}
#
# example: maker_run.py 2 10
########################################################################################################################
import os
import sys
import shutil
import inspect

# sys.path.append(r'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Acan_NEFF/NEFF_maker/MakerScripts')
import makerFunctions as Funcs
import makerSettings as Sets


if __name__ == '__main__':

    pass_num = int(sys.argv[1])
    threads = int(sys.argv[2])

    x = [1]

    for i in x:
        new_paths = Sets.Paths(i=i, pass_num=pass_num)
        old_paths = Sets.Paths(i=i, pass_num=(pass_num - 1))
        my_exes = Sets.Executables()

        # Makes MAKER directory and cd's to it
        Funcs.make_dir(new_paths.maker_dir)
        os.chdir(new_paths.maker_dir)

        # Create the 3 MAKER2 .ctl files
        my_exes.maker(arg='-CTL')

        with open('mod_maker_exe.ctl', 'w') as new_exe_file, open('maker_exe.ctl', 'r') as exe_file:
            for line in exe_file:
                line = line.strip()
                if line.split('=')[0] == 'RepeatMasker':
                    line = line.split('=')[0] + '=' + my_exes.rep_mask_exe + ' #' + line.split('#')[1]
                new_exe_file.write(line + '\n')

        # Edit maker_opts.ctl - This file contatins the options for the maker run. Add more elif's if you want to add
        # more options
        with open(f'maker_opts_edit.ctl', 'w') as new_opt_file, open('maker_opts.ctl', 'r') as opt_file:
            for line in opt_file:
                line = line.strip()
                line_begin = line.split('=')[0]
                try:
                    line_end = line.split('#')[1]
                except IndexError:
                    pass

                if line_begin == 'genome':
                    line = f'{line_begin}={old_paths.genome_path} #{line_end}'
                elif line_begin == 'est_pass' and pass_num < 4:
                    line = f'{line_begin}=1 #{line_end}'
                elif line_begin == 'est2genome' and pass_num == 1:
                    line = f'{line_begin}=1 #{line_end}'
                # elif line_begin == 'protein2genome' and pass_num == 1:
                #     line = f'{line_begin}=1 #{line_end}'
                elif line_begin == 'keep_preds' and pass_num == 1:
                    line = f'{line_begin}=1 #{line_end}'
                elif line_begin == 'single_exon' and pass_num == 1:
                    line = f'{line_begin}=1 #{line_end}'
                # elif line_begin == 'est' and pass_num < 4:
                #     line = f'{line_begin}={old_paths.rna_seq} #{line_end}'
                elif line_begin == 'altest' and pass_num < 4:
                    line = f'{line_begin}={old_paths.rna_seq} #{line_end}'
                elif line_begin == 'rmlib' and pass_num < 4:
                    line = f'{line_begin}={old_paths.rep_mod_out} #{line_end}'
                elif line_begin == 'gmhmm' and pass_num < 4:
                    line = f'{line_begin}={old_paths.gm_hmm} #{line_end}'
                elif line_begin == 'maker_gff' and 1 < pass_num < 4:
                    line = f'{line_begin}={old_paths.maker_all_gff} #{line_end}'
                elif line_begin == 'snaphmm' and 1 < pass_num < 4:
                    line = f'{line_begin}={old_paths.snap_hmm} #{line_end}'
                elif line_begin == 'augustus_species' and 1 < pass_num < 4:
                    line = f'{line_begin}={old_paths.base_name} #{line_end}'
                elif line_begin == 'pred_gff' and pass_num == 4:
                    line = f'{line_begin}={old_paths.maker_all_gff} #{line_end}'
                elif line_begin == 'keep_preds' and pass_num == 4:
                    line = f'{line_begin}=1 #{line_end}'
                elif line_begin == 'model_org' and pass_num == 4:
                    line = f'{line_begin}= #{line_end}'
                elif line_begin == 'rmlib' and pass_num == 4:
                    line = f'{line_begin}= #{line_end}'
                elif line_begin == 'repeat_protein' and pass_num == 4:
                    line = f'{line_begin}= #{line_end}'

                new_opt_file.write(line + '\n')

        # Removes orginal, unedited maker_opts.ctl and maker_exe.ctl
        os.remove('maker_opts.ctl')
        os.remove('maker_exe.ctl')

        snap_line = ''
        busco_line = ''
        next_pass = ''

        if pass_num < 3:
            # Commands to train SNAP
            # these lines go in the job script that follows if pass 1 or 2
            Funcs.make_dir(new_paths.snap_dir)
            snap_line = f'''        
                cd {new_paths.snap_dir}
                maker2zff {new_paths.maker_all_gff}
                fathom -categorize 1000 genome.ann genome.dna
                fathom -export 1000 -plus uni.ann uni.dna 
                forge export.ann export.dna
                hmm-assembler.pl {new_paths.base_name}_snap . > {new_paths.base_name}_snap.hmm
            '''

            # BUSCO run after passes 1 and 2
            # Trains AUGUSTUS
            Funcs.make_dir(new_paths.busco_dir)
            busco_line = f'''
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

            # Line to to contiue to next pass
            next_pass = f'python3 ~/Scripts/MakerScripts/maker_run.py {pass_num+1} {threads}'

        if pass_num == 3:
            # Final BUSCO run (likely needs to change to pass_num == 3)
            # fourth pass through maker likely not necesarry
            # Does NOT train AUGUSTUS
            Funcs.make_dir(new_paths.final_busco_dir)
            busco_line = f'''
                cd {new_paths.final_busco_dir}
                export AUGUSTUS_CONFIG_PATH={new_paths.aug_config_path}
                python {my_exes.busco_exe} \\
                    -i {new_paths.maker_prot} \\
                    -l {new_paths.busco_lineage} \\
                    -o {new_paths.base_name}_prot \\
                    -m proteins \\
                    -c {threads}
            '''

        # Job Script
        # Runs maker followed by BUSCO and Snap
        # only trains Snap after passes 1 and 2
        # BUSCO only trains AUGUSTUS after passes 1 and 2
        with open(f'maker.sh', 'w') as job_script:
            job_str = inspect.cleandoc(f'''\
                #!/bin/bash
                #$ -S /bin/bash
                . /etc/profile
                #$ -cwd
                #$ -pe threads {threads}

                conda activate maker
                module load maker/2.31
                
                cd ~/Scripts/MakerScripts
                python3 -c "import makerFunctions; makerFunctions.get_aug_species({i}, {pass_num-1})"
                
                cd {new_paths.maker_dir}
                maker -c {threads - 1} -fix_nucleotides maker_opts_edit.ctl maker_bopts.ctl mod_maker_exe.ctl

                cd {new_paths.maker_out}
                fasta_merge -d {new_paths.maker_log}
                gff3_merge -d {new_paths.maker_log}

                {snap_line}
                
                {busco_line}
                                
                ''')
            job_script.write(job_str)
        # Submit job script
        my_exes.qsub(arg=f'maker.sh')


main()

