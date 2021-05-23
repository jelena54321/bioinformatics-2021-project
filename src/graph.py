from Bio import SeqIO, Seq, SeqRecord
from node import Node, Overlap
import pafpy
import random
import copy


SEQ_ID_MIN = 0.9
LEN_DELTA = 1000
NUM_ELEMENTS = 3
NUM_RUNS = 1
DELTA_TRIALS = 100
MAX_SEQ_LEN = 2_000_000
MAX_SEQ_BETWEEN_LEN = 150_000
MIN_LEN = 10_000


class SortHelper:

    def __init__(self, node, score):
        self.node = node
        self.score = score


class Graph:

    COMPARATOR_BY_SCORE_AND_LENGTH = lambda el: (el.score, el.node.len)
    OVERLAP_SCORE = lambda *args: Graph.overlap_score(args[0])

    def __init__(self, contigs, reads):
        self.contigs = contigs
        self.reads = reads

    @staticmethod
    def construct(reads_to_contigs, reads_to_reads):
        contigs = dict()
        reads = dict()

        Graph.process_overlaps(
            reads_to_contigs, contigs, reads,
            Graph.should_add_overlap
        )
        print_str = (
            'graph.Graph.construct >> '
            'Generated overlaps between contigs and reads.'
        )
        print(print_str)

        Graph.process_overlaps(
            reads_to_reads, reads, reads,
            lambda *_: True
        )
        print('graph.Graph.construct >> Generated overlaps between reads.')

        return Graph(contigs, reads)

    def generate_paths(self):
        all_paths = []
        all_path_lens = []
        for run in range(NUM_RUNS):
            print_str = (
                'graph.Graph.generate_paths >> '
                f'Finding Paths - run: {run + 1}/{NUM_RUNS}'
            )
            print(print_str)

            contigs = list(self.contigs.values())
            first_node = random.choice(contigs)
            while len(first_node.nodes) == 0:
                first_node = random.choice(contigs)

            for node in first_node.nodes:
                paths, path_lens = self.dfs(
                    node, [first_node],
                    Graph.OVERLAP_SCORE,
                    self.next_using_the_best_score
                )
                all_paths.extend(paths)
                all_path_lens.extend(path_lens)

                paths, path_lens = self.dfs(
                    node, [first_node],
                    Graph.extension_score,
                    self.next_using_the_best_score
                )
                all_paths.extend(paths)
                all_path_lens.extend(path_lens)

            n_next_nodes = len(first_node.nodes)
            for i in range(n_next_nodes + DELTA_TRIALS):
                paths, path_lens = self.dfs(
                    first_node, [],
                    Graph.extension_score,
                    self.next_using_monte_carlo
                )
                all_paths.extend(paths)
                all_path_lens.extend(path_lens)

        print(f'graph.Graph.generate_paths >> Found {len(all_paths)} paths.')
        return all_paths, all_path_lens

    def generate_sequence(self, paths, path_lens, contigs, reads, out):
        biggest_group = max(Graph.generate_groups(paths, path_lens), key=len)
        best_path = Graph.best_path(biggest_group)

        used_nodes = set(map(lambda node: node.id, best_path))

        Graph.load_sequences(contigs, self.contigs, used_nodes)
        print('graph.Graph.generate_sequence >> Loaded contigs.')

        Graph.load_sequences(reads, self.reads, used_nodes)
        print('graph.Graph.generate_sequence >> Loaded reads.')

        seq = Graph.build_sequence(best_path)
        record = SeqRecord.SeqRecord(seq, id='output_sequence')
        with open(out, 'w') as handle:
            SeqIO.write([record], handle, 'fasta')

        print_str = (
            'graph.Graph.generate_sequence >> '
            'Written generated sequence in the output file.'
        )
        print(print_str)

    @staticmethod
    def process_overlaps(path, queries, targets, should_add_overlap):

        with pafpy.PafFile(path) as paf:
            for ol in paf:
                if ol.is_unmapped(): continue
                if Graph.sequence_identity(ol) < SEQ_ID_MIN: continue

                if ol.qname not in queries:
                    q = Node(ol.qname, ol.qlen)
                    queries[q.id] = q
                    compl_q = Node.complement(q)
                    queries[compl_q.id] = compl_q
                else:
                    q = queries[ol.qname]
                    compl_q = queries[Node.complement_id(q.id)]

                if ol.tname not in targets:
                    t = Node(ol.tname, ol.tlen)
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

                if not should_add_overlap(ol, t, compl_t, queries, targets):
                    continue

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

    def get_node_complement(self, node):
        nodes = self.reads if node.id in self.reads else self.contigs
        return nodes[Node.complement_id(node.id)]

    def dfs(self, start_node, start_path, scores_fn, next_fn):
        paths = []
        path_lens = []

        open_nodes_between = [[start_node]]
        open_paths = [copy.copy(start_path)]
        start_path_len = 0 if len(start_path) == 0 else start_path[0].len
        open_path_lens = [start_path_len]

        contigs = set(self.contigs.values())
        visited = set()

        while len(open_nodes_between) != 0:

            node = open_nodes_between[0].pop()
            active_path = open_paths[0]
            
            #ako mi nije counting
            ext_len = 0

            if node not in contigs:
                if len(active_path) == 0:
                    ext_len = node.len
                else:
                    ol = active_path[-1].overlap_for_node(node)
                    ext_len = node.len - ol.tend

            path_len = open_path_lens[0] + ext_len

            if len(open_nodes_between[0]) == 0:
                open_nodes_between.pop()
                open_paths.pop()
                open_path_lens.pop()

            if node in visited:
                continue
            else:
                visited.add(node)

            node_compl = self.get_node_complement(node)
            if node in active_path or node_compl in active_path: continue

            if len(node.nodes) == 0 or path_len > MAX_SEQ_BETWEEN_LEN: continue

            is_read = node not in contigs
            connected_contig = set(node.nodes).intersection(contigs)
            if is_read and len(connected_contig) != 0:
                contig = connected_contig.pop()

                path = copy.copy(active_path)
                path.append(node)
                path.append(contig)
                paths.append(path)

                ol = node.overlap_for_node(contig)
                ext_len = contig.len - ol.tend
                path_lens.append(path_len + ext_len)

                open_nodes_between.clear()
                open_paths.clear()
                open_path_lens.clear()

                open_nodes_between.append(next_nodes)
                open_paths = [copy.copy(next_nodes)]
                open_path_lens = [next_nodes]

                continue

            next_nodes = self.next_nodes(node, active_path, scores_fn, next_fn)
            if len(next_nodes) == 0: continue

            path = copy.copy(active_path)
            path.append(node)

            open_nodes_between.append(next_nodes)
            open_paths.append(path)
            open_path_lens.append(path_len)

        return paths, path_lens

    def next_nodes(self, node, path, scores_fn, next_fn):
        scores = []
        for i in range(len(node.nodes)):
            next_node = node.nodes[i]
            ol = node.overlaps[i]

            score = scores_fn(ol, node, next_node)
            scores.append(SortHelper(next_node, score))

        return next_fn(scores, path)

    def next_using_the_best_score(self, scores, path):
        scores = sorted(
            scores,
            key=Graph.COMPARATOR_BY_SCORE_AND_LENGTH,
            reverse=True
        )

        next = []

        # TODO: CHECK
        # if this is the best approach
        for i in range(NUM_ELEMENTS):
            node = scores[i].node
            if node in path or self.get_node_complement(node) in path:
                continue

            next.append(node)

        return next

    def next_using_monte_carlo(self, scores, path):
        next = []
        nodes = []
        weights = []
        for node in scores:
            if node.node in path or self.get_node_complement(node.node) in path:
                continue

            nodes.append(node.node)
            weights.append(node.score)

        if len(nodes) == 0: return next

        node = random.choices(nodes, weights)
        next.append(node[0])
        return next

    @staticmethod
    def generate_groups(paths, path_lens):
        if all(path_len < MIN_LEN for path_len in path_lens):
            return [paths]

        len_groups = dict()

        for i in range(len(paths)):
            path = paths[i]
            path_len = path_lens[i]

            found_group = False
            for len_group in len_groups:
                if ((path_len >= len_group - LEN_DELTA) and
                    (path_len <= len_group + LEN_DELTA)):

                    len_groups[len_group].append(path)
                    found_group = True
                    break

            if not found_group:
                len_groups[path_len] = [path]

        return list(len_groups.values())

    @staticmethod
    def best_path(group):
        max_ol_score = -1
        best_path = None
        for path in group:
            ol_score = 0
            for i in range(len(path) - 1):
                ol = path[i].overlap_for_node(path[i+1])
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
            ol = node.overlap_for_node(path[i+1])

            seq += node.seq[start:ol.qstart]
            start = ol.tstart

        seq += path[-1].seq[start:]
        return seq

    @staticmethod
    def load_sequences(path, nodes, used_nodes):
        with open(path) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if record.id in used_nodes:
                    nodes[record.id].seq = record.seq
                    continue

                compl_id = Node.complement_id(record.id)
                if compl_id in used_nodes:
                    nodes[compl_id].seq = Node.complement_sequence(record.seq)

    @staticmethod
    def sequence_identity(ol: pafpy.PafRecord):
        return ol.mlen / ol.blen

    @staticmethod
    def average_overlap_length(ol):
        return (ol.qend-ol.qstart+ol.tend-ol.tstart) / 2

    @staticmethod
    def overlap_score(ol):
        if type(ol) is pafpy.PafRecord:
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
        return ol_score + ext_len/2 - (oh_1 + oh_2)/2

    @staticmethod
    def should_add_overlap(ol, read, compl_read, contigs, reads):
        if len(read.nodes) == 0 and len(compl_read.nodes) == 0: return True

        r = read if len(read.nodes) > 0 else compl_read
        old_contig = r.nodes[0]
        old_overlap = r.overlaps[0]

        old_ol_score = Graph.overlap_score(old_overlap)
        ol_score = Graph.overlap_score(ol)
        if old_ol_score > ol_score: return False

        r.remove_overlap(old_contig)

        compl_r = reads[Node.complement_id(r.id)]
        compl_c = contigs[Node.complement_id(old_contig.id)]
        compl_c.remove_overlap(compl_r)

        return True
