from typing import List
import pafpy
from Bio import SeqIO, Seq, SeqRecord
from node import Node, Overlap
import random

SEQ_ID_MIN = 0.85
LEN_DELTA = 100
HEURISTIC_DEPTH = 30
NUM_ELEMENT = 3
NUM_RUNS = 10
NUM_TRIALS = 3

class Graph:

    def __init__(self, contigs, reads):
        self.contigs = contigs
        self.reads = reads

    @staticmethod
    def construct(reads_to_contigs, reads_to_reads):
        contigs = dict()
        reads = dict()

        Graph.process_overlaps(reads_to_contigs, contigs, reads, Graph.should_add_overlap)
        print('graph.Graph.construct >> Generated overlaps between contigs and reads.')

        Graph.process_overlaps(reads_to_reads, reads, reads, lambda *_: True)
        print('graph.Graph.construct >> Generated overlaps between reads.')

        return Graph(contigs, reads)


    #klasa koja ce mi pomoc sortirati
    class helper_sort:
        score: int
        node: Node
        def __init__(self, n:Node, s:int):
            self.score = s
            self.node = n


    #stavlja nejbolje elemente na next
    def get_best_paths(self, l, next, current_path):
        for i in l:
            if i.node in current_path:
                continue
            next.append(i.node)

            if len(next) == NUM_ELEMENT:
                break

    def monte_carlo(self, all_extension_scores, next, current_path):
        nodes = []
        weight = []
        for node in all_extension_scores:
            if node in current_path: continue

            nodes.append(node.node)
            weight.append(node.score)

        if len(nodes) == 0: return

        node = random.choices(nodes, weight)[0]
        next.append(node)

        return

    def dfs(self, current_node: Node, heuristic_depth: int, all_found_paths, current_path, funID):
        #found new counting
        if current_node in self.contigs:
            heuristic_depth = HEURISTIC_DEPTH

        #dfs dosao do kraj, ili je zbog heuristike ili nea vise djece
        if heuristic_depth == 0 or len(current_node.nodes) == 0:
            all_found_paths.append(current_path)
            return

        # dobi sve elemente
        all_overlap_scores = []
        all_extension_scores = []
        for i in range(len(current_node.nodes)):
            next_node = current_node.nodes[i]
            ol = current_node.overlaps[i]

            ol_score = Graph.overlap_score(ol)
            all_overlap_scores.append(helper_sort(next_node, ol_score))

            ext_score = Graph.extension_score(ol, current_node, next_node)
            all_extension_scores.append(helper_sort(next_node, ext_score))

        all_overlap_scores = sorted(all_overlap_scores, key=lambda el: el.score, reverse=True)
        all_extension_scores = sorted(all_extension_scores, key=lambda el: el.score, reverse=True)

        next = []

        if funID == 0:
            self.get_best_paths(all_overlap_scores, next, current_path)
        elif funID == 1:
            self.get_best_paths(all_extension_scores, next, current_path)
        elif funID == 2:
            self.monte_carlo(all_extension_scores, next, current_path)

        for i in next:
            current_path.append(i)
            self.dfs(i, heuristic_depth - 1, all_found_paths, current_path, funID)
            #backtracking
            current_path.pop()


    def generate_paths(self):
        """
        TODO: Traverse through graph using depth first search and return list of found
        paths (list of node lists).
        """
        current_path = []
        all_found_paths = []
        for run in range(NUM_RUNS):
            print(f'graph.Graph.generate_paths >> Finding Paths: {run + 1}')

            first_node = random.choice(list(self.contigs.values()))
            while(len(first_node.nodes) == 0):
                first_node = random.choice(list(self.contigs.values()))

            current_path.append(first_node)
            for node in first_node.nodes:
                current_path.append(node)
                self.dfs(node, HEURISTIC_DEPTH, all_found_paths, current_path, 0)
                self.dfs(node, HEURISTIC_DEPTH, all_found_paths, current_path, 1)
                for i in range(NUM_TRIALS):
                    self.dfs(node, HEURISTIC_DEPTH, all_found_paths, current_path, 2)
                current_path.pop()

        return all_found_paths


    def generate_sequence(self, paths, contigs, reads, out):
        biggest_group = max(Graph.generate_groups(paths), key=len)
        best_path = best_path(biggest_group)

        Graph.load_sequences(contigs, self.contigs)
        print('graph.Graph.generate_sequence >> Loaded contigs.')

        Graph.load_sequences(reads, self.reads)
        print('graph.Graph.generate_sequence >> Loaded reads.')

        seq = Graph.build_sequence(best_path)
        record = SeqRecord.SeqRecord(Seq.Seq(seq), id='output_sequence')
        with open(out_path, 'w') as handle:
            SeqIO.write([record], handle, 'fasta')

        print('graph.Graph.generate_sequence >> Written generated sequence in the output file.')

    @staticmethod
    def process_overlaps(path, queries, targets, should_add_overlap):

        with pafpy.PafFile(path) as paf:
            for ol in paf:
                if ol.is_unmapped(): continue
                if Graph.sequence_identity(ol) < SEQ_ID_MIN: continue

                if ol.qname not in queries:
                    q = Node(ol.qname)
                    queries[q.id] = q
                    compl_q = Node.complement(q)
                    queries[compl_q.id] = compl_q
                else:
                    q = queries[ol.qname]
                    compl_q = queries[Node.complement_id(q.id)]

                if ol.tname not in targets:
                    t = Node(ol.tname)
                    targets[t.id] = t
                    compl_t = Node.complement(t)
                    targets[compl_t.id] = compl_t
                else:
                    t = targets[ol.tname]
                    compl_t = targets[Node.complement_id(t.id)]

                q_left = ol.qstart
                q_right = ol.qlen - ol.qend
                t_left = ol.tstart
                t_right = ol.tlen - ol.tend

                q_contains_t = t_left <= q_left and t_right <= q_right
                t_contains_q = q_left <= t_left and q_right <= t_right
                if q_contains_t or t_contains_q: continue

                if not should_add_overlap(ol, t, compl_t, queries, targets): continue

                seq_id = Graph.sequence_identity(ol)
                if ol.strand == pafpy.Strand.Forward:
                    if t_left > q_left:
                        overlap = Overlap(
                            ol.tstart, ol.tend,
                            ol.qstart, ol.qend,
                            seq_id
                        )
                        t.add_overlap(q, overlap)

                        overlap = Overlap(
                            ol.qlen - ol.qend, ol.qlen - ol.qstart,
                            ol.tlen - ol.tend, ol.tlen - ol.tstart,
                            seq_id
                        )
                        compl_q.add_overlap(compl_t, overlap)

                    else:
                        overlap = Overlap(
                            ol.qstart, ol.qend,
                            ol.tstart, ol.tend,
                            seq_id
                        )
                        q.add_overlap(t, overlap)

                        overlap = Overlap(
                            ol.tlen - ol.tend, ol.tlen - ol.tstart,
                            ol.qlen - ol.qend, ol.qlen - ol.qstart,
                            seq_id
                        )
                        compl_t.add_overlap(compl_q, overlap)

                else:
                    if t_left > q_right:
                        overlap = Overlap(
                            ol.qstart, ol.qend,
                            ol.tlen - ol.tend, ol.tlen - ol.tstart,
                            seq_id
                        )
                        q.add_overlap(compl_t, overlap)

                        overlap = Overlap(
                            ol.tstart, ol.tend,
                            ol.qlen - ol.qend, ol.qlen - ol.qstart,
                            seq_id
                        )
                        t.add_overlap(compl_q, overlap)

                    else:
                        overlap = Overlap(
                            ol.tlen - ol.tend, ol.tlen - ol.tstart,
                            ol.qstart, ol.qend,
                            seq_id
                        )
                        compl_t.add_overlap(q, overlap)

                        overlap = Overlap(
                            ol.qlen - ol.qend, ol.qlen - ol.qstart,
                            ol.tstart, ol.tend,
                            seq_id
                        )
                        compl_q.add_overlap(t, overlap)

    @staticmethod
    def generate_groups(paths):
        len_groups = dict()

        for path in paths:
            path_len = len(path)

            found_group = False
            for len_group in len_groups:
                if (path_len >= len_group - LEN_DELTA) and (path_len <= len_group + LEN_DELTA):
                    len_groups[len_group].append(path)
                    found_group = True
                    break

            if not found_group:
                len_groups[path_len] = [path]

        return len_groups.values()

    @staticmethod
    def best_path(group):
        max_ol_score = -1
        best_path = None
        for path in group:
            ol_score = 0
            for i in range(len(path) - 1):
                ol = path[i].overlap_for_node(path[i + 1])
                ol_score += Graph.overlap_score(ol)

            if ol_score > max_ol_score:
                max_ol_score = ol_score
                best_path = path

        return best_path

    @staticmethod
    def build_sequence(path):
        seq = ''
        start = 0
        for i in range(len(path) - 1):
            node = path[i]
            ol = node.overlap_for_node(path[i + 1])

            seq += node.seq[start:ol.qstart]
            start = ol.tstart

        seq += path[-1].seq[start:]

        return seq

    @staticmethod
    def load_sequences(path, nodes):
        with open(path) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                nodes[record.id].seq = record.seq
                nodes[Node.complement_id(record.id)].seq = Node.complement_sequence(record.seq)

    @staticmethod
    def sequence_identity(ol: pafpy.PafRecord):
        return ol.mlen / ol.blen

    @staticmethod
    def average_overlap_length(ol):
        return (ol.qend - ol.qstart + ol.tend - ol.tstart) / 2

    @staticmethod
    def overlap_score(ol):
        if isinstance(ol, pafpy.PafRecord):
            seq_id = Graph.sequence_identity(ol)
        else:
            seq_id = ol.seq_id
        return seq_id * Graph.average_overlap_length(ol)

    @staticmethod
    def extension_score(ol: Overlap, q, t):
        ol_score = Graph.overlap_score(ol)
        oh_1 = len(q) - ol.qend
        oh_2 = ol.tstart
        ext_len = len(t) - ol.tend
        return ol_score + ext_len / 2 - (oh_1 + oh_2) / 2

    @staticmethod
    def should_add_overlap(ol, read, compl_read, contigs, reads):
        if len(read.nodes) == 0 and len(compl_read.nodes) == 0: return True

        r = read if len(read.nodes) > 0 else compl_read
        old_contig = r.nodes[0]
        old_overlap = r.overlaps[0]

        if Graph.overlap_score(old_overlap) > Graph.overlap_score(ol): return False

        r.remove_overlap(old_contig)

        compl_r = reads[Node.complement_id(r.id)]
        compl_c = contigs[Node.complement_id(old_contig.id)]
        compl_c.remove_overlap(compl_r)

        return True
