import sys
from Bio import AlignIO
from Bio.AlignIO.MafIO import MafWriter
from collections import defaultdict
#from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
def startPosition(sequence):
        if sequence.annotations["strand"] == -1:
                return sequence.annotations["srcSize"] - sequence.annotations["start"] - sequence.annotations["size"]
        return sequence.annotations["start"]

def addBlocks(infile, outfile, sequences):
        genomes = defaultdict(list)
        genome_size = defaultdict(int)
        old_records, new_records = [], []
        for i, mafblock in enumerate(AlignIO.parse(infile, "maf")):
                old_records.append(mafblock)
                for sequence in mafblock:
                        genome_size[sequence.id] = sequence.annotations["srcSize"]
                        start = startPosition(sequence)
                        end = start + sequence.annotations["size"]
                        genomes[sequence.id].append((start, end))
        for genome, path in genomes.items():
                path.sort(key = lambda t: t[0])
                if path[0][0] > 0:
                        record = SeqRecord(Seq(sequences[genome][0:path[0][0]]), id=genome, annotations = {"srcSize": genome_size[genome], "start": 0, "strand": 1, "size": path[0][0]})
                        aln = MultipleSeqAlignment([record])
                        new_records.append(aln)
                for j in range(len(path)-1):
                        end = path[j][1]
                        start = path[j+1][0]
                        if end < start:
                                record = SeqRecord(Seq(sequences[genome][end:start]), id=genome, annotations = {"srcSize": genome_size[genome], "start": end, "strand": 1, "size": start-end})
                                aln = MultipleSeqAlignment([record])
                                new_records.append(aln)
                if path[-1][1] < genome_size[genome]:
                        record = SeqRecord(Seq(sequences[genome][path[-1][1]:]), id=genome, annotations = {"srcSize": genome_size[genome], "start": path[-1][1], "strand": 1, "size": genome_size[genome]-path[-1][1]})
                        aln = MultipleSeqAlignment([record])
                        new_records.append(aln)
        handle = open(outfile, "w")
        mw = MafWriter(handle)
        mw.write_header()
        new_records += old_records
        for multiple_alignment in new_records:
                mw.write_alignment(multiple_alignment)
        handle.close()

def main():
        infile, outfile, seq_file = sys.argv[1], sys.argv[2], sys.argv[3]
        sequences = defaultdict(str)
        name = None
        with open(seq_file) as f:
                for line in f:
                        line = line.strip()
                        if line.startswith(">"):
                                name = line[1:]
                        else:
                                sequences[name] = line
        addBlocks(infile, outfile, sequences)

if __name__ == "__main__":
        main()
