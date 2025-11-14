class Walk:
        def __init__(self, start, end, orient, length):
                self.start = start
                self.end = end
                self.orient = orient
                self.length = length

class Genome:
        def __init__(self, idx, start):
                self.id = idx
                self.start = start
                self.end = 0
                self.name = ""
                self.size = 0

class Graph:
        def __init__(self, in_file):
                self.vtx_seq = ["",]
                self.paths = [0,]
                self.genomes = []
                self.used = [False,]

                g_id, pos = 0,1
                with open(in_file) as file:
                        for line in file:
                                if line.startswith("S"):
                                        line = line.strip().split()[2]
                                        self.vtx_seq.append(line.upper())
                                elif line.startswith("P"):
                                        genome = Genome(g_id, pos)
                                        line = line.strip().split()
                                        genome.name = line[1]
                                        path = line[2].split(',')
                                        for el in path:
                                                vtx = int(el[:-1])
                                                strand = 1 if el[-1]=='+' else -1
                                                self.paths.append(vtx*strand)
                                                genome.size += len(self.vtx_seq[vtx])
                                                pos += 1
                                        genome.end = pos-1
                                        pos += 1
                                        self.genomes.append(genome)
                                        self.paths.append(0)
                                        g_id += 1
                self.occurrences = [[] for i in range(len(self.vtx_seq))]
                self.used = [False for i in range(len(self.paths))]
                self.set_occurrences()

        def set_occurrences(self):
                for i, el in enumerate(self.paths):
                        if el!=0:
                                self.occurrences[abs(el)].append(i)

        def find_seeds(self, vtx_id, PARAM_a):
                collinear_seeds = []
                carrying_seed_orientation = None
                if len(self.occurrences[vtx_id])<=PARAM_a:
                        for occ in self.occurrences[vtx_id]:
                                if self.used[occ]==False:
                                        if not collinear_seeds:
                                                orient = 1
                                                carrying_seed_orientation = 1 if self.paths[occ]>0 else -1
                                        else:
                                                orient_on_path = 1 if self.paths[occ]>0 else -1
                                                orient = orient_on_path*carrying_seed_orientation
                                        collinear_seeds.append(Walk(occ, occ, orient, len(self.vtx_seq[vtx_id])))
                return collinear_seeds, carrying_seed_orientation
