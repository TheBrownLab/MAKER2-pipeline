########################################################################################################################
# Robert E. Jones (rej110@msstate.edu)
# GitHub Repo: https://github.com/rej110/myMAKER2-Pipeline
# usage:   maker_run.py -i PATH_TO_INFILE
#          run remove_ambiguity.py -h for options
# Details: Removes IUPAC amiguity codes from nucleotide sequences (in FASTA format) and replaces them with N.
#          Some programs, such as RepeatModeler, cannot handle ambiguous nucleotides.
########################################################################################################################
import argparse
from Bio import SeqIO
from Bio.Seq import MutableSeq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Removes IUPAC ambiguity codes and replaces them with N',
                                     usage="remove_ambiguity.py -i fasta_file")
    parser.add_argument('-i', '--in_file', type=str, required=True,
                        help='path to file in fasta format')
    args = parser.parse_args()

    standard = ['A', 'a', 'T', 't', 'C', 'c', 'G', 'g']

    with open(args.in_file, 'r') as infile, open('test.fas', 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            my_seq = MutableSeq(str(record.seq))

            for i in range(0, len(my_seq)):
                if my_seq[i] not in standard:
                    my_seq[i] = 'N'

            outfile.write(f'>{record.description}\n')
            outfile.write(f'{str(my_seq)}\n')
