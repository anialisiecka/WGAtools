from collections import defaultdict

def sort_vtxs(graph):
        vtx_ids_order = defaultdict(int)
        for vtx in range(1,len(graph.vtx_seq)):
                if len(graph.occurrences[vtx])>1:
                        vtx_ids_order[vtx] = len(graph.occurrences[vtx])*len(graph.vtx_seq[vtx])
        vtx_ids_order = sorted(vtx_ids_order.keys(), key=lambda x: vtx_ids_order[x], reverse=True)
        return vtx_ids_order

def walk_length(start, end, graph):
        w_len = 0
        pos = start
        while (pos<=end):
                w_len += len(graph.vtx_seq[abs(graph.paths[pos])])
                pos += 1
        return w_len

def mark_vtxs(graph, block):
        for walk in block.collinear_walks:
                for i in range(walk.start, walk.end+1):
                        graph.used[i] = True

def unmark_vtxs(graph, block):
        for walk in block.collinear_walks:
                for i in range(walk.start, walk.end+1):
                        graph.used[i] = False

def scoring_function(block, graph, PARAM_b):
        score = 0
        for walk_id, walk in enumerate(block.collinear_walks):
                if block.scores[walk_id].q1>PARAM_b or block.scores[walk_id].q3>PARAM_b:
                        return -1
                score += walk.length - (block.scores[walk_id].q1 + block.scores[walk_id].q3)
        return score

def remove_short_walks(block, PARAM_m):
        new_walks, new_match_carrying, new_scores = [], [], []
        for i, walk in enumerate(block.collinear_walks):
                if walk.length >= PARAM_m:
                        new_walks.append(walk)
                        new_match_carrying.append(block.match_carrying[i])
                        new_scores.append(block.scores[i])
        block.collinear_walks = new_walks
        block.match_carrying = new_match_carrying
        block.scores = new_scores

complement_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
        seq_complement = ""
        for s in seq[::-1]:
                seq_complement += complement_dict[s]
        return seq_complement
