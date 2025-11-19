import sys
from Bio import AlignIO
from collections import defaultdict

def main():
        input = sys.argv[1]
        coverage = 0.0
        genomes_len = defaultdict(int)
        for mafblock in AlignIO.parse(input, "maf"):
                if len(mafblock)>1:
                        for sequence in mafblock:
                                genomes_len[sequence.id] = sequence.annotations["srcSize"]
                                coverage += sequence.annotations["size"]
        total_len = 0.0
        for key, val in genomes_len.items():
                total_len += val
        print("coverage", coverage/total_len)

if __name__ == "__main__":
        main()
