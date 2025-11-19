import sys
from Bio import AlignIO
from collections import defaultdict

def main():
        input_maf = sys.argv[1]
        average_identity = 0.0
        col_numbers = 0
        for mafblock in AlignIO.parse(input_maf, "maf"):
                n = len(mafblock)
                if n>1:
                        columns = [defaultdict(int) for i in range(len(mafblock[0].seq))]
                        for sequence in mafblock:
                                for i,c in enumerate(sequence.seq):
                                        if c != '-':
                                                columns[i][c.upper()] += 1
                        #col_numbers += len(mafblock[0].seq)
                        for d in columns:
                                if d:
                                        col_numbers += 1
                                        for key, val in d.items():
                                                average_identity += val*(val-1)/(n*(n-1))
        print("average identity score", average_identity/col_numbers)

if __name__ == "__main__":
        main()
