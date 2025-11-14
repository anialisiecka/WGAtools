import sys
import argparse
from graph import Graph
from block_tools import sort_vtxs
from block_detector import find_block

sys.path.append('poapy')
from poa import poa_align
def walk_start(walk, graph):
        start=0
        genome_id = None
        for i, genome in enumerate(graph.genomes):
                if genome.start <= walk.start <= genome.end:
                        genome_id=i
                        if walk.orient == 1:
                                if walk.start==genome.start:
                                        return i, 0
                                else:
                                        for pos in range(genome.start, walk.start):
                                                start += len(graph.vtx_seq[abs(graph.paths[pos])])
                                        return i, start
                        else:
                                for pos in range(genome.start, walk.end+1):
                                        start += len(graph.vtx_seq[abs(graph.paths[pos])])
                                start = genome.size-start
                                return i, start
        return

def wga(in_file, out_file, _match, _mismatch, _gap, a, b, m):
        graph = Graph(in_file)
        vtx_ids_order = sort_vtxs(graph)
        count = 0
        maf_file = open(out_file, 'w')
        maf_file.write('##maf version=1 scoring=tba.v8\n\n')
        for vtx_id in vtx_ids_order:
                block = find_block(graph, vtx_id, PARAM_a=a, PARAM_b=b, PARAM_m=m)
                if block is not None:
                        maf_file.write("a\n")
                        alignment = poa_align(block, graph, _match=_match, _mismatch=_mismatch, _gap=_gap)[1:]
                        for i, walk in enumerate(block.collinear_walks):
                                genome_id, start = walk_start(walk, graph)
                                w_orient = "+" if walk.orient==1 else "-"
                                maf_file.write("s\t"+graph.genomes[genome_id].name+"\t"+str(start)+"\t"+str(walk.length)+"\t"+w_orient+"\t"+str(graph.genomes[genome_id].size)+"\t"+alignment[i][1]+"\n")
                        maf_file.write("\n")
        maf_file.close()

if __name__=="__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', required=True, help='Path to input .gfa file.')
        parser.add_argument('-o', required=True, help='Path to output .maf file.')
        parser.add_argument('-a', default=150, required=False, help='Abundance pruning parameter, default to 150. Vertices occurring more than a times are not considered as collinear walk seeds.')
        parser.add_argument('-b', default=200, required=False, help='Maximal bubble size, default to 200 residues.')
        parser.add_argument('-m', default=50, required=False, help='Minimal collinear walk length, default to 50.')
        parser.add_argument('--match', default=5, required=False, help='Match score in alignment. Default to 5.')
        parser.add_argument('--mismatch', default=-4, required=False, help='Mismatch penalty in alignment. Default to -4.')
        parser.add_argument('--gap', default=-8, required=False, help='Gap penalty in alignment. Default to -8.')
        args = parser.parse_args()

        wga(args.i, args.o, args.match, args.mismatch, args.gap, a=int(args.a), b=int(args.b), m=int(args.m))
