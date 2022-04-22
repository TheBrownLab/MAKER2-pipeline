########################################################################################################################
# Robert E. Jones (rej110@msstate.edu)
# GitHub Repo: https://github.com/rej110/myMAKER2-Pipeline
# usage:   summarize.py [options]
#          run summarize.py -h for options
# Details: This is a script that summarizes the results of maker_run.py
########################################################################################################################
import re
import os


def get_basename():
    cwd = os.getcwd()
    _base = cwd.split('/')[-1]

    return _base


def get_busco_stats():
    busco = []
    for i in range(1, 5):
        busco_file = f'./BUSCO_Pass{i}/run_{base}.genome_prot/short_summary_{base}.genome_prot.txt'

        with open(busco_file, 'r') as infile:
            for line in infile:
                line = line.strip()
                complete = re.findall(r'C:[0-9]{0,2}\.[0-9]', line)
                if complete:
                    busco.append(float(complete[0].split(':')[1]))

    return ",".join(repr(e) for e in busco)


def get_gene_num():
    genes = []
    for i in range(1, 5):
        genes_file = f'./Maker_Pass{i}/{base}.genome.maker.output/{base}.genome.all.maker.proteins.fasta'
        gene_count = 0

        with open(genes_file, 'r') as infile:
            for line in infile:
                line = line.strip()
                gene_count += len(re.findall(r'>', line))
            genes.append(gene_count)

    return ",".join(repr(e) for e in genes)


if __name__ == '__main__':
    base = get_basename()

    with open(f'{base}.maker.summary.csv', 'w') as outfile:
        outfile.write(f'BUSCO Scores,{get_busco_stats()}\n'
                      f'Gene Count,{get_gene_num()}\n')

