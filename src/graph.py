from Bio import SeqIO, Seq, SeqRecord
from node import Node, Overlap
import pafpy
import random
import copy


SEQ_ID_MIN = 0.9
LEN_DELTA = 100_000
NUM_ELEMENTS = 2
DELTA_TRIALS = 300
MAX_BETWEEN_LEN = 250_000
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
            'Added overlaps between contigs and reads.'
        )
        print(print_str)

        Graph.process_overlaps(
            reads_to_reads, reads, reads,
            lambda *_: True
        )
        print('graph.Graph.construct >> Added overlaps between reads.')

        return Graph(contigs, reads)

    def generate_paths(self):
        paths = []
        path_lens = []

        print('graph.Graph.generate_paths >> Finding paths...')

        for contig in list(self.contigs.values()):

            n_successors = len(contig.nodes)
            if n_successors == 0: continue

            for node in contig.nodes:
                path, path_len = self.dfs(
                    node, [contig],
                    Graph.OVERLAP_SCORE,
                    self.next_using_the_best_score
                )
                if path is None: continue
                paths.append(path)
                path_lens.append(path_len)

                path, path_len = self.dfs(
                    node, [contig],
                    Graph.extension_score,
                    self.next_using_the_best_score
                )
                if path is None: continue
                paths.append(path)
                path_lens.append(path_len)

            for _ in range(n_successors + DELTA_TRIALS):
                path, path_len = self.dfs(
                    contig, [],
                    Graph.extension_score,
                    self.next_using_monte_carlo
                )
                if path is None: continue
                paths.append(path)
                path_lens.append(path_len)

        print(f'graph.Graph.generate_paths >> Found {len(paths)} paths.')
        return paths, path_lens

    def generate_sequence(self, paths, path_lens, contigs, reads, out):
        paths, path_lens = self.filter_unique_paths(paths, path_lens)
        paths, path_lens = self.remove_subpaths(paths, path_lens)

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
        cpt_path = None
        cpt_path_len = 0

        open_nodes = [[start_node]]
        open_paths = [copy.copy(start_path)]
        start_path_len = 0 if len(start_path) == 0 else start_path[0].len
        open_path_lens = [start_path_len]
        open_between_lens = [0]

        contigs = set(self.contigs.values())
        visited = set()
        if len(start_path) != 0:
            visited.add(start_path[0])

        while len(open_nodes) != 0:

            node = open_nodes[0].pop(0)
            active_path = open_paths[0]

            between_len = open_between_lens[0]
            is_read = node not in contigs
            if len(active_path) == 0:
                ext_len = node.len
            else:
                prev_node = active_path[-1]
                ol = prev_node.overlap_for_node(node)
                ext_len = -(prev_node.len-ol.qstart) + (node.len-ol.tstart)

                if is_read:
                    between_len += node.len - ol.tend

            path_len = open_path_lens[0] + ext_len

            if len(open_nodes[0]) == 0:
                open_nodes.pop(0)
                open_paths.pop(0)
                open_path_lens.pop(0)
                open_between_lens.pop(0)

            if node in visited:
                continue
            else:
                visited.add(node)

            node_compl = self.get_node_complement(node)
            if node in active_path or node_compl in active_path: continue

            if len(node.nodes) == 0 or between_len > MAX_BETWEEN_LEN:
                continue

            connected_contig = set(node.nodes).intersection(contigs)
            if is_read and len(connected_contig) != 0:
                contig = connected_contig.pop()
                contig_compl = self.get_node_complement(contig)
                if contig in active_path or contig_compl in active_path:
                    continue

                cpt_path = copy.copy(active_path)
                cpt_path.append(node)
                tmp_path = copy.copy(cpt_path)
                cpt_path.append(contig)

                ol = node.overlap_for_node(contig)
                ext_len = -(node.len-ol.qstart) + (contig.len-ol.tstart)
                cpt_path_len = path_len + ext_len

                open_nodes.clear()
                open_paths.clear()
                open_path_lens.clear()
                open_between_lens.clear()
                visited.clear()

                open_nodes.insert(0, [contig])
                open_paths.insert(0, tmp_path)
                open_path_lens.insert(0, path_len)
                open_between_lens.insert(0, 0)

                continue

            next_nodes = self.next_nodes(node, active_path, scores_fn, next_fn)
            if len(next_nodes) == 0: continue

            path = copy.copy(active_path)
            path.append(node)

            open_nodes.insert(0, next_nodes)
            open_paths.insert(0, path)
            open_path_lens.insert(0, path_len)
            open_between_lens.insert(0, between_len)

        return cpt_path, cpt_path_len

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
        for i in range(min(NUM_ELEMENTS, len(scores))):
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

    def filter_unique_paths(self, paths, path_lens):
        to_be_filtered = set()

        n_paths = len(paths)
        for i in range(n_paths):
            path = paths[i]

            for j in range(i + 1, n_paths):
                other_path = paths[j]

                path_len = len(path)
                if path_len != len(other_path): continue

                if path[0] == other_path[0]:
                    is_reversed = False
                elif path[0] == self.get_node_complement(other_path[-1]):
                    is_reversed = True
                else:
                    continue

                different = False
                for k in range(1, path_len):
                    node = path[k]

                    if is_reversed:
                        idx = path_len - k - 1
                        other_node = self.get_node_complement(other_path[idx])
                    else:
                        other_node = other_path[k]

                    if node != other_node:
                        different = True
                        break

                if not different:
                    to_be_filtered.add(j)

        unique_paths = []
        lens = []
        for i in range(n_paths):
            if i in to_be_filtered: continue

            unique_paths.append(paths[i])
            lens.append(path_lens[i])

        print_str = (
            'graph.Graph.filter_unique_paths >> '
            f'Removed {n_paths - len(unique_paths)} path duplicates.'
        )
        print(print_str)

        return unique_paths, lens

    def remove_subpaths(self, paths, path_lens):
        to_be_filtered = set()

        n_paths = len(paths)
        for i in range(n_paths):
            path = paths[i]
            path_ctgs = self.filter_contigs(path)

            for j in range(i + 1, n_paths):
                other_path = paths[j]
                other_path_ctgs = self.filter_contigs(other_path)

                if len(path_ctgs) == len(other_path_ctgs): continue

                path_ctgs_sublist_of_other_path_ctgs = Graph.is_sublist(
                    path_ctgs,
                    other_path_ctgs
                )
                other_path_ctgs_sublist_of_path_ctgs = Graph.is_sublist(
                    other_path_ctgs,
                    path_ctgs
                )

                if path_ctgs_sublist_of_other_path_ctgs:
                    to_be_filtered.add(i)
                elif other_path_ctgs_sublist_of_path_ctgs:
                    to_be_filtered.add(j)

        unique_paths = []
        lens = []
        for i in range(n_paths):
            if i in to_be_filtered: continue

            unique_paths.append(paths[i])
            lens.append(path_lens[i])

        print_str = (
            'graph.Graph.remove_subpaths >> '
            f'Removed {n_paths - len(unique_paths)} subpaths.'
        )
        print(print_str)

        return unique_paths, lens

    def filter_contigs(self, path):
        contigs = set(self.contigs.values())
        return list(filter(lambda n: n in contigs, path))

    @staticmethod
    def is_sublist(ls1, ls2):
        if len(ls1) > len(ls2): return False

        for i in range(len(ls2)):
            if ls1[0] != ls2[i]: continue

            j = 1
            while j < len(ls1) and i+j < len(ls2) and ls1[j] == ls2[i+j]:
                j += 1

            if j == len(ls1): return True

        return False

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

            ol_score /= len(path) - 1

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
            for record in SeqIO.parse(handle, Graph.determine_type(path)):
                if record.id in used_nodes:
                    nodes[record.id].seq = record.seq
                    continue

                compl_id = Node.complement_id(record.id)
                if compl_id in used_nodes:
                    nodes[compl_id].seq = Node.complement_sequence(record.seq)

    @staticmethod
    def determine_type(path):
        return 'fastq' if path.endswith('.fastq') else 'fasta'

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
