import re
import sys

sys.path.append(r'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/MakerPipeline/Acan_NEFF/NEFF_maker/MakerScripts')
import makerSettings as Sets


def main():
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

            outfile.write(f'ref{i},{Sets.key[i]},{complete_list[0]},{complete_list[1]},{complete_list[2]}\n')


main()
