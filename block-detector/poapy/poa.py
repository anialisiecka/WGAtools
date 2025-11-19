import poagraph
import seqgraphalignment
from collections import namedtuple

complement_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
AlignmentTuple = namedtuple('AlignmentTuple', ['sequence', 'stringidxs', 'nodeidxs'])

def reverse_complement(seq):
        seq_complement = ""
        for s in seq[::-1]:
                seq_complement += complement_dict[s]
        return seq_complement

def find_edges_between(graph, start, end, walks_to_poa):
        node = graph.nodedict[start]
        walks_to_align = {w_id for e in node.outEdges for w_id in node.outEdges[e].labels}
        node = graph.nodedict[end]
        walks_to_align = walks_to_align.intersection({w_id for e in node.inEdges for w_id in node.inEdges[e].labels})
        to_use = set(range(start+1, end))
        for w_id in walks_to_align:
                found_start = False
                if w_id == -1:
                        continue
                for nid in walks_to_poa[w_id]:
                        if nid==end:
                                break
                        if nid==start:
                                found_start=True
                        elif found_start:
                                to_use.add(nid)
        return to_use

def carrying_to_poa_idx(block, vg):
        carrying_to_poa = {}
        carrying_seq = ""
        poa_v_nr = 0
        for i, v_oriented in enumerate(block.carrying_path):
                carrying_seq += vg.vtx_seq[v_oriented] if v_oriented>0 else reverse_complement(vg.vtx_seq[-v_oriented])
                carrying_to_poa[i] = [poa_v_nr]
                poa_v_nr = len(carrying_seq)
                carrying_to_poa[i].append(poa_v_nr-1)
        return carrying_to_poa, carrying_seq

def align_matching_fragment(start_nid, end_nid, sequence, seq_len_so_far, old_alignment):
        if old_alignment is None:
                old_alignment = [sequence, list(range(len(sequence))), list(range(start_nid, end_nid+1))]
        else:
                old_alignment[0] += sequence
                old_alignment[1] += [sidx+seq_len_so_far for sidx in range(len(sequence))]
                old_alignment[2] += list(range(start_nid, end_nid+1))
        return old_alignment

def add_to_alignment(old_alignment, new_alignment, seq_len_so_far):
        old_alignment[0] += new_alignment.sequence
        old_alignment[1] += [None if sidx is None else sidx+seq_len_so_far for sidx in new_alignment.stringidxs]
        old_alignment[2] += new_alignment.nodeidxs
        return old_alignment

def poa_align(block, vg, _match, _mismatch, _gap, globalAlign=True, simple=True):
        carrying_to_poa, carrying_seq = carrying_to_poa_idx(block, vg)
        graph = poagraph.POAGraph(carrying_seq, -1) # initialize graph with carrying_path
        walks_to_poa = {}
        entire_sequences = []

        for w_id, walk in enumerate(block.collinear_walks):
                matches = block.match_carrying[w_id]
                to_start = walk.start if walk.orient==1 else walk.end
                to_end = walk.end+1 if walk.orient==1 else walk.start-1

                if len(matches)==1:
                        match_carrying = matches[0][0]
                        start_nid, end_nid = carrying_to_poa[match_carrying]
                        entire_sequence = ""
                        for pos in range(to_start, to_end, walk.orient):
                                v_oriented = vg.paths[pos]*walk.orient
                                entire_sequence += vg.vtx_seq[v_oriented] if v_oriented>0 else reverse_complement(vg.vtx_seq[-v_oriented])
                        old_alignment = align_matching_fragment(start_nid, end_nid, entire_sequence, 0, None)
                        entire_sequences.append(entire_sequence)
                        old_alignment = AlignmentTuple(*old_alignment)
                        walks_to_poa[w_id] = graph.incorporateSeqAlignment(old_alignment, entire_sequence, label=w_id)
                        continue

                first = True
                old_alignment = None
                i_on_match = 0
                sequence, entire_sequence = "", ""
                pos = to_start
                while pos != to_end:
                        v_oriented = vg.paths[pos]*walk.orient
                        if matches[i_on_match][1] != pos: # there is no match on this position of walk
                                sequence += vg.vtx_seq[v_oriented] if v_oriented>0 else reverse_complement(vg.vtx_seq[-v_oriented])
                                pos += walk.orient
                        elif sequence == "":
                                match_carrying = matches[i_on_match][0]
                                start_nid, end_nid = carrying_to_poa[match_carrying]
                                sequence = vg.vtx_seq[v_oriented] if v_oriented>0 else reverse_complement(vg.vtx_seq[-v_oriented])
                                old_alignment = align_matching_fragment(start_nid, end_nid, sequence, len(entire_sequence), old_alignment)
                                i_on_match += 1
                                pos += walk.orient
                                entire_sequence += sequence
                                sequence = ""
                                first = False
                        else:
                                last_match = matches[i_on_match-1]
                                match = matches[i_on_match]
                                start_nid = carrying_to_poa[last_match[0]][1]
                                end_nid = carrying_to_poa[match[0]][0]
                                to_use = find_edges_between(graph, start_nid, end_nid, walks_to_poa)
                                if to_use:
                                        alignment = seqgraphalignment.SeqGraphAlignment(sequence, graph, fastMethod=not simple, globalAlign=globalAlign, matchscore=_match, mismatchscore=_mismatch, gapsco>
                                else:
                                        alignment = AlignmentTuple(sequence, range(len(sequence)), [None for _ in range(len(sequence))])
                                old_alignment = add_to_alignment(old_alignment, alignment, len(entire_sequence))
                                entire_sequence += sequence
                                sequence = ""
                entire_sequences.append(entire_sequence)
                old_alignment = AlignmentTuple(*old_alignment)
                walks_to_poa[w_id] = graph.incorporateSeqAlignment(old_alignment, entire_sequence, label=w_id)
        alignments = graph.generateAlignmentStrings()
        return alignments
