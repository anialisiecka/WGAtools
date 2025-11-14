from Bio import AlignIO
from collections import defaultdict
from graph import Graph
import sys

def intToStrand(s):
        if s == 1:
                return "+"
        return "-"

def start_position(sequence):
        if sequence.annotations["strand"] == -1:
                return sequence.annotations["srcSize"] - sequence.annotations["start"] - sequence.annotations["size"]
        return sequence.annotations["start"]

def firstNonGappedColumn(sequence):
        for i, residue in enumerate(sequence):
                if residue != "-":
                        return i, residue.upper()

def lastNonGappedColumn(sequence):
        for i in range(len(sequence)-1, -1, -1):
                residue = sequence[i]
                if residue != "-":
                        return i, residue.upper()

def poa(genome, sequence, strand, block_id, graph):
        if strand == 1:
                prev_col, prev_residue = firstNonGappedColumn(sequence)
                prev_node = (block_id, prev_col, prev_residue)
                prev_node_id = graph.update_nodes(prev_node)
                j = prev_col + 1
                for i in range(j, len(sequence)):
                        if sequence[i] != "-":
                                next_node = (block_id, i, sequence[i].upper())
                                next_node_id = graph.update_nodes(next_node)
                                edge = (prev_node_id, strand, next_node_id, strand)
                                graph.update_edges(edge, edge, genome)
                                prev_node_id = next_node_id
        else:
                prev_col, prev_residue = lastNonGappedColumn(sequence)
                prev_node = (block_id, prev_col, prev_residue)
                prev_node_id = graph.update_nodes(prev_node)
                j = prev_col - 1
                for i in range(j, -1, -1):
                        if sequence[i] != "-":
                                next_node = (block_id, i, sequence[i].upper())
                                next_node_id = graph.update_nodes(next_node)
                                genomic_edge = (prev_node_id, strand, next_node_id, strand)
                                edge = genomic_edge[2:]+ genomic_edge[:2]
                                graph.update_edges(genomic_edge, edge, genome)
                                prev_node_id = next_node_id

def wga2vg(infile, outfile):
        paths = defaultdict(list)
        G = Graph()
        for i, mafblock in enumerate(AlignIO.parse(infile, "maf")):
                for sequence in mafblock:
                        seq_id = sequence.id
                        paths[seq_id].append((i, start_position(sequence), sequence))
        for seq_id, path in paths.items():
                path.sort(key = lambda t: t[1])
                prev_block = path[0]
                prev_strand = prev_block[2].annotations["strand"]
                prev_seq = prev_block[2].seq
                if prev_strand == 1:
                        prev_col, prev_residue = lastNonGappedColumn(prev_seq)
                else:
                        prev_col, prev_residue = firstNonGappedColumn(prev_seq)
                prev_block_id = prev_block[0]
                prev_node = (prev_block_id, prev_col, prev_residue)
                prev_node_id = G.update_nodes(prev_node)
                poa(seq_id, prev_seq, prev_strand, prev_block_id, G)
                for j in range(1, len(path)):
                        next_block = path[j]
                        next_block_id = next_block[0]
                        next_strand = next_block[2].annotations["strand"]
                        next_seq = next_block[2].seq
                        if next_strand == 1:
                                next_col, next_residue = firstNonGappedColumn(next_seq)
                        else:
                                next_col, next_residue = lastNonGappedColumn(next_seq)
                        next_node = (next_block_id, next_col, next_residue)
                        next_node_id = G.update_nodes(next_node)
                        genomic_edge = (prev_node_id, prev_strand, next_node_id, next_strand)
                        if prev_block_id < next_block_id:
                                edge = (prev_node_id, prev_strand, next_node_id, next_strand)
                        elif prev_block_id > next_block_id:
                                edge = (next_node_id, next_strand, prev_node_id, prev_strand)
                        else:
                                if prev_col < next_col:
                                        edge = (prev_node_id, prev_strand, next_node_id, next_strand)
                                else:
                                        edge = (next_node_id, next_strand, prev_node_id, prev_strand)
                        G.update_edges(genomic_edge, edge, seq_id)
                        poa(seq_id, next_seq, next_strand, next_block_id, G)
                        prev_seq = next_seq
                        prev_block_id = next_block_id
                        prev_node_id = next_node_id
                        prev_col = next_col
                        prev_residue = next_residue
                        prev_strand = next_strand
                        prev_node = next_node
        out = open(outfile, "w")
        out.write('H'+'\t'+'VN:Z:1.0'+'\n')
        for node, node_id in G.nodes.items():
                residue = node[2]
                out.write("S" + "\t" + str(node_id) + "\t" + residue + "\n")
        for edge in G.edges:
                out.write("L" + "\t" + str(edge[0]) + "\t" + intToStrand(edge[1]) + "\t" + str(edge[2]) + "\t" + intToStrand(edge[3]) + "\t" + "0M" + "\n")
        for seq_id, path in G.genomic_paths.items():
                path_to_string = ""
                edge = path[0]
                path_to_string = str(edge[0]) + intToStrand(edge[1]) + "," + str(edge[2]) + intToStrand(edge[3])
                for j in range(1, len(path)):
                        edge = path[j]
                        path_to_string += "," + str(edge[2]) + intToStrand(edge[3])
                out.write("P" + "\t" + seq_id + "\t" + path_to_string + "\t" + "*" + "\n")
        out.close()

def main():
        infile, outfile = sys.argv[1], sys.argv[2]
        wga2vg(infile, outfile)

if __name__ == "__main__":
        main()
