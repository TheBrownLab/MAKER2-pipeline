from Bio import SeqIO


if __name__ == '__main__':
    key = []
    with open('/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Genomes/RenamedGenomesKey2.csv') as infile:
        for line in infile:
            line = line.strip()
            key.append(line.split(','))

    for item in key:
        genome = f'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Genomes/{item[2]}'
        pseudo_contig = f'/mnt/scratch/brownlab/rej110/AcanthamoebaGenomes/Maker_alt_EST/ref{item[0]}/data/ref{item[0]}.genome.fas'

        with open(genome, 'r') as ref, open(pseudo_contig, 'w') as outfile:
            records = SeqIO.parse(ref, 'fasta')
            new_contig = ''
            for record in records:
                list_of_strings = record.format('fasta').split('\n')[1:]
                new_contig += (('N'*600) + "".join(list_of_strings).strip())

            outfile.write(f'>PsuedoContig_with_600N_bewtween_contigs\n{new_contig}')

