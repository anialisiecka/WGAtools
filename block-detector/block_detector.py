import copy
from graph import Walk
from block import Block, Extensions, Score
from block_tools import scoring_function, mark_vtxs, unmark_vtxs, remove_short_walks

def reverse_block(block):
        block.carrying_path.reverse()
        for i in range(len(block.carrying_path)):
                block.carrying_path[i] *= (-1)
        num_of_vtxs = len(block.carrying_path)-1
        for i in range(len(block.match_carrying)):
                matches = reversed(block.match_carrying[i])
                block.match_carrying[i] = [(-m[0]+num_of_vtxs,m[1]) for m in matches]
        for walk in block.collinear_walks:
                walk.orient *= (-1)

def minus(collinear_seed):
        seeds = []
        for s in collinear_seed:
                seeds.append(Walk(s.start, s.end, -s.orient, s.length))
        return seeds

def find_block(graph, vtx_id, PARAM_a, PARAM_b, PARAM_m):
        collinear_seeds, carrying_seed_orientation = graph.find_seeds(vtx_id, PARAM_a)
        if not collinear_seeds:
                return None
        new_seeds = minus(collinear_seeds)
        new_block = Block(vtx_id, len(graph.vtx_seq[vtx_id]), new_seeds, -carrying_seed_orientation)
        mark_vtxs(graph, new_block)
        new_score = scoring_function(new_block, graph, PARAM_b)
        best_score = new_score
        best_block = copy.deepcopy(new_block)
        while new_score >= 0:
                w0 = new_block.carrying_path[-1] # oriented last vtx on carrying path
                Q = Extensions(new_block.collinear_walks, graph, w0, PARAM_b)
                r = Q.get_carrying_path_extension(graph, new_block)
                if not r:
                        break
                for w in r:
                        new_block.carrying_path.append(w)
                        update_walks(new_block, w, Q, graph, PARAM_b)
                        new_score = scoring_function(new_block, graph, PARAM_b)
                        new_block.carrying_path_length += len(graph.vtx_seq[abs(w)])
                        if new_score < 0:
                                break
                        if new_score > best_score:
                                best_block = copy.deepcopy(new_block)
                                best_score = new_score
        reverse_block(best_block)
        unmark_vtxs(graph, new_block)
        new_block = copy.deepcopy(best_block)
        new_score = best_score
        mark_vtxs(graph, best_block)
        while new_score >= 0:
                # find carrying path extension
                w0 = new_block.carrying_path[-1] # oriented last vtx on carrying path
                Q = Extensions(new_block.collinear_walks, graph, w0, PARAM_b)
                r = Q.get_carrying_path_extension(graph, new_block)
                if not r:
                        break
                for w in r:
                        new_block.carrying_path.append(w)
                        update_walks(new_block, w, Q, graph, PARAM_b)
                        new_score = scoring_function(new_block, graph, PARAM_b)
                        new_block.carrying_path_length += len(graph.vtx_seq[abs(w)])
                        if new_score < 0:
                                break
                        if new_score > best_score:
                                best_block = copy.deepcopy(new_block)
                                best_score = new_score
        if best_score > 0:
                remove_short_walks(best_block, PARAM_m)
                best_block.remove_double_matches()
                if len(best_block.collinear_walks)>1:
                        unmark_vtxs(graph, new_block)
                        mark_vtxs(graph, best_block)
                        return best_block
                unmark_vtxs(graph, new_block)
        return None

def find_walk_to_extend(block, extensions, graph, occ, w):
        walk_to_extend = None
        for i, extension in extensions.extensions.items():
                walk = block.collinear_walks[i]
                if extension[0]<=occ<=extension[1] and walk.orient*graph.paths[occ]==w:
                        if walk_to_extend is None:
                                walk_to_extend = i
                        else:
                                if walk.orient==1 and walk.end>block.collinear_walks[walk_to_extend].end:
                                        walk_to_extend = i
                                elif walk.orient==-1 and walk.start<block.collinear_walks[walk_to_extend].start:
                                        walk_to_extend = i
        return walk_to_extend

def update_walks(block, w, extensions, graph, PARAM_b):
        updated_scores = set([])
        for occ in graph.occurrences[abs(w)]:
                if graph.used[occ]==False: # if occurrence is not used
                        walk_to_extend = find_walk_to_extend(block, extensions, graph, occ, w)
                        if walk_to_extend is not None:
                                walk = block.collinear_walks[walk_to_extend]
                                if walk.orient==1:
                                        for i in range(walk.end+1, occ+1):
                                                graph.used[i] = True
                                                walk.length += len(graph.vtx_seq[abs(graph.paths[i])])
                                        walk.end = occ
                                        block.scores[walk_to_extend].q3 = 0
                                else:
                                        for i in range(occ, walk.start):
                                                graph.used[i] = True
                                                walk.length += len(graph.vtx_seq[abs(graph.paths[i])])
                                        walk.start = occ
                                        block.scores[walk_to_extend].q1 = 0
                                updated_scores.add(walk_to_extend)
                                block.match_carrying[walk_to_extend].append((len(block.carrying_path)-1, occ))
                                extensions.update_extension(walk, graph, walk_to_extend, PARAM_b)
                        else:
                                orient_on_carrying_path = 1 if w>0 else -1
                                orient_on_path = 1 if graph.paths[occ]>0 else -1
                                orient = orient_on_carrying_path*orient_on_path
                                graph.used[occ] = True
                                block.collinear_walks.append(Walk(occ, occ, orient, len(graph.vtx_seq[abs(w)])))
                                if orient==1:
                                        block.scores.append(Score(block.carrying_path_length, 0))
                                else:
                                        block.scores.append(Score(0, block.carrying_path_length))
                                updated_scores.add(len(block.collinear_walks)-1)
                                extensions.update_extension(block.collinear_walks[-1], graph, len(block.collinear_walks)-1, PARAM_b)
                                block.match_carrying.append([(len(block.carrying_path)-1, occ)])
        for w_id, walk in enumerate(block.collinear_walks):
                if w_id not in updated_scores:
                        if walk.orient==1:
                                block.scores[w_id].q3 += len(graph.vtx_seq[abs(w)])
                        else:
                                block.scores[w_id].q1 += len(graph.vtx_seq[abs(w)])
