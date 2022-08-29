from Bio import SeqIO
import argparse as ap


if __name__ == '__main__':
    
    parser = ap.ArgumentParser()
    parser.add_argument(type=str, dest='in_file', help='Path to input genome')
    parser.add_argument(type=str, dest='out_file', help='Path to output psudo-contig')
    args = parser.parse_args()

    genome = args.in_file
    pseudo_contig = args.out_file

    with open(genome, 'r') as ref, open(pseudo_contig, 'w') as outfile:
        records = SeqIO.parse(ref, 'fasta')
        new_contig = ''
        for record in records:
            list_of_strings = record.format('fasta').split('\n')[1:]
            new_contig += (('N'*600) + "".join(list_of_strings).strip())

        outfile.write(f'>PsuedoContig_with_600N_bewtween_contigs\n{new_contig}')

