from collections import defaultdict

class Score:
        def __init__(self, q1, q3):
                self.q1 = q1
                self.q3 = q3
class Path:
        def __init__(self, start, end, distance, walk_id):
                self.start = start
                self.end = end
                self.distance = distance
                self.walk_id = walk_id

class Block:
        def __init__(self, seed_idx, seed_length, seed_occurrences, carrying_seed_orientation):
                self.carrying_path = [seed_idx*carrying_seed_orientation]
                self.collinear_walks = seed_occurrences
                self.scores = [Score(0,0) for w in self.collinear_walks]
                self.carrying_path_length = seed_length
                self.match_carrying = [[(0, w.start)] for w in self.collinear_walks]

        def remove_double_matches(self):
                for w_id in range(len(self.collinear_walks)):
                        matches = self.match_carrying[w_id].copy()
                        for i in range(len(matches)-1, 0, -1):
                                prev_match, _ = self.match_carrying[w_id][i-1]
                                match, _ = self.match_carrying[w_id][i]
                                if prev_match==match:
                                        if i==1:
                                                self.match_carrying[w_id].pop(i)
                                        else:
                                                self.match_carrying[w_id].pop(i-1)

class Extensions:
        def __init__(self, collinear_walks, graph, w0, PARAM_b):
                self.extensions = {}
                self.coverage = defaultdict(int)
                self.shortest_walk = {}

                for walk_id, walk in enumerate(collinear_walks):
                        if walk.orient==1:
                                pos = walk.end
                                proximal = walk.end+1
                        else:
                                pos = walk.start
                                proximal = walk.start-1
                        if graph.used[proximal] or graph.paths[proximal]==0:
                                continue
                        p_length = 0
                        i = proximal
                        u = graph.paths[pos]*walk.orient # last vtx of a walk
                        while True:
                                if graph.used[i] or graph.paths[i]==0:
                                        i -= walk.orient
                                        break
                                v_oriented = graph.paths[i]*walk.orient
                                self.coverage[v_oriented] += 1
                                p_length += len(graph.vtx_seq[abs(v_oriented)])
                                if u==w0:
                                        if v_oriented not in self.shortest_walk or self.shortest_walk[v_oriented].distance > p_length:
                                                self.shortest_walk[v_oriented] = Path(pos, i, p_length, walk_id)
                                if p_length >= PARAM_b:
                                        break
                                i += walk.orient
                        self.extensions[walk_id] = (min(proximal,i), max(proximal,i))

        def get_carrying_path_extension(self, graph, block):
                highest_coverage = -1
                w0_to_t = None
                for v_oriented in self.shortest_walk:
                        if self.coverage[v_oriented] > highest_coverage:
                                w0_to_t = self.shortest_walk[v_oriented] #Path(start, end, length, walk_id)
                                highest_coverage = self.coverage[v_oriented]

                if w0_to_t is None:
                        return []

                r = []
                walk = block.collinear_walks[w0_to_t.walk_id]
                for i in range(w0_to_t.start+walk.orient, w0_to_t.end+walk.orient, walk.orient):
                        r.append(graph.paths[i]*walk.orient)
                return r

        def update_extension(self, walk, graph, walk_id, PARAM_b):
                if walk.orient == 1:
                        proximal = walk.end + 1
                else:
                        proximal = walk.start -1
                if graph.paths[proximal]==0 or graph.used[proximal]==True:
                        if walk_id in self.extensions:
                                del self.extensions[walk_id]
                else:
                        distal = proximal
                        p_length = 0
                        while True:
                                if graph.used[distal]==True or graph.paths[distal]==0:
                                        distal -= walk.orient
                                        break
                                p_length += len(graph.vtx_seq[abs(graph.paths[distal])])
                                if p_length >= PARAM_b:
                                        break
                                distal += walk.orient
                        self.extensions[walk_id] = (min(proximal,distal), max(proximal,distal))
