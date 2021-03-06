import random
import copy
import pafpy
from Bio import SeqIO, SeqRecord
from graph_models import Node, Overlap
from search_models import SortHelper, SearchState, TestData


SEQ_ID_MIN = 0.9
LEN_DELTA = 100_000
NUM_ELEMENTS = 2
DELTA_TRIALS = 300
MAX_GAP_LEN = 250_000
MIN_LEN = 10_000


class Graph:
    '''Class that represents an overlap graph.

    :param contigs: Contig nodes
    :param reads: Read nodes
    '''

    COMPARATOR_BY_SCORE_AND_LENGTH = lambda el: (el.score, el.node.len)
    OVERLAP_SCORE = lambda *args: Graph.overlap_score(args[0])

    def __init__(self, contigs, reads):
        self.contigs = contigs
        self.reads = reads

    @staticmethod
    def construct(reads_to_contigs, reads_to_reads):
        '''Builds a graph based on provided PAF align files.

        :param reads_to_contigs: PAF file with reads aligned to contigs
        :param reads_to_reads: PAF file with reads aligned to other reads
        :returns: New Graph instance
        '''

        contigs = dict()
        reads = dict()

        Graph.process_overlaps(
            reads_to_contigs, contigs, reads,
            Graph.should_add_overlap
        )
        print_str = (
            '\ngraph.Graph.construct >> '
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
        '''Generates node paths as a result of a depth first search.

        All contigs are used as path starting points. In order to avoid
        traversing through whole graph, this method uses three different
        heuristics:
        (1) Expanding current node with overlapping node that has the
        highest overlap score
        (2) Expanding current node with overlapping node that has the
        highest extension score
        (3) Expanding current node with randomly chosen overlapping node.
        The probability of a connecting node being selcted is proportional
        to its extension score.

        :returns: List of found paths and their corresponding lengths
        '''

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
        '''Generates final sequence based on the list of found paths.

        Method filters only unique paths and discards subpaths. Remaining
        paths are grouped based on their lengths. From the biggest group
        the path with the highest overlap score is chosen as a final path.
        Using data from the provided files, result sequence is built and
        written to the output file.

        :param paths: List of found paths
        :param path_lens: List of corresponding path lengths
        :param contigs: FASTA/FASTQ file with contig data
        :param reads: FASTA/FASTQ file with reads data
        :param out: FASTA file where final sequence will be written
        '''

        if len(paths) == 0:
            print('graph.Graph.generate_sequence >> No paths. Exiting.')
            return None

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

        return self.prepare_test_data(best_path)

    @staticmethod
    def process_overlaps(path, queries, targets, should_add_overlap):
        '''Processes PAF align files to build a graph.

        Method builds overlaps between two reads or reads and contigs
        based on the data found in the PAF file.

        Unmapped aligns, aligns with sequence identity smaller than the
        threshold value and aligns with reads that are contained in
        other sequences are discarded.

        :param path: PAF file with align records
        :param queries: Dictionary where record queries will be saved
        :param targets: Dictionary where record targets will be saved
        :param should_add_overlap: Function that determines whether
        overlap will be added to graph
        '''

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
        '''Gets node complement in a graph.

        :param node: Node
        :returns: Node complement
        '''

        nodes = self.reads if node.id in self.reads else self.contigs
        return nodes[Node.complement_id(node.id)]

    def dfs(self, start_node, start_path, scores_fn, next_fn):
        '''Runs depth first search through the graph.

        Search does not allow paths where distance between two contigs
        is bigger than the threshold value. When method encounters
        a contig, search continues as though said contig is the starting
        point.

        :param start_node: Node that will be visited first
        :param start_path: Initial path
        :param scores_fn: Function for determining overlap scores
        :param next_fn: Function for determining node successors
        :returns: Found node path
        '''

        cpt_path = None
        cpt_path_len = 0

        start_state = SearchState(
            [start_node],
            copy.copy(start_path),
            0 if len(start_path) == 0 else start_path[0].len,
            0
        )
        open = [start_state]

        contigs = set(self.contigs.values())
        visited = set()
        if len(start_path) != 0:
            visited.add(start_path[0])

        while len(open) != 0:

            state = open[-1]

            node = state.nodes.pop()
            active_path = state.path
            gap_len = state.gap_len

            is_read = node not in contigs
            if len(active_path) == 0:
                ext_len = node.len
            else:
                prev_node = active_path[-1]
                ol = prev_node.overlap_for_node(node)
                ext_len = -(prev_node.len-ol.qstart) + (node.len-ol.tstart)

                if is_read:
                    gap_len += node.len - ol.tend

                    if len(active_path) > 1:
                        gap_len -= (prev_node.len-ol.qend)

            path_len = state.path_len + ext_len

            if len(state.nodes) == 0:
                open.pop()

            if node in visited: continue

            visited.add(node)

            node_compl = self.get_node_complement(node)
            if node in active_path or node_compl in active_path: continue

            if len(node.nodes) == 0 or gap_len > MAX_GAP_LEN:
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

                open.clear()
                open.append(SearchState([contig], tmp_path, path_len, 0))

                continue

            next_nodes = self.next_nodes(node, active_path, scores_fn, next_fn)
            if len(next_nodes) == 0: continue

            path = copy.copy(active_path)
            path.append(node)

            open.append(SearchState(next_nodes[::-1], path, path_len, gap_len))

        return cpt_path, cpt_path_len

    def next_nodes(self, node, path, scores_fn, next_fn):
        '''Gets children nodes based on the scores_fn and next_fn.

        :param node: Node that will be expanded
        :param path: Active node path
        :param scores_fn: Function for determining node scores
        :param next_fn: Function for determining next nodes
        :returns: Next nodes to be visited
        '''

        scores = []
        for i in range(len(node.nodes)):
            next_node = node.nodes[i]
            ol = node.overlaps[i]

            score = scores_fn(ol, node, next_node)
            scores.append(SortHelper(next_node, score))

        return next_fn(scores, path)

    def next_using_the_best_score(self, scores, path):
        '''Gets at most NUM_ELEMENTS elements with the best score.

        :param scores: List of SortHelper instances
        :param path: Active node path
        :return: Next nodes to be visited using the best score
        '''

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
        '''Gets next node using Monte Carlo method.

        Method randomly chooses next node. The probability of a node
        being chosen is proportional to the node extension score.

        :param scores: List of SortHelper instances
        :param path: Active node path
        :returns: Next node to be visited using Monte Carlo
        '''

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
        '''Filters unique paths.

        :param paths: List of paths
        :param path_lens: List of corresponding path lengths
        :returns: Unique paths
        '''

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
        '''Removes subpaths.

        Subpaths are all paths that contain contigs that are found in
        other paths. Order of contigs is taken into account.

        :param paths: List of paths
        :param path_lens: List of corresponding path lengths
        :returns: Paths without subpaths
        '''

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
        '''Returns list of contigs found in the provided path.

        :param path: Path
        :returns: List of contigs included in the path
        '''

        contigs = set(self.contigs.values())
        return list(filter(lambda n: n in contigs, path))

    @staticmethod
    def is_sublist(ls1, ls2):
        '''Determines whether one is sublist of another.

        :param ls1: First list
        :param ls2: Second list
        :returns: True if ls1 is sublist of ls2
        '''

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
        '''Groups paths into groups.

        Paths with simmilar lengths are grouped into the same group.
        If all path lengths are smaller than MIN_LEN, then all paths
        are put in the same group.

        :param path: Paths list
        :param path_lens: Path lengths list
        :returns: Grouped paths
        '''

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
        '''Returns best path according to the highest overlap score

        In order to prevent paths being chosen only based on their
        lengths, total overlap score for a path is devided by the
        number of present overlaps.

        :param group: Group of paths
        :returns: Path with the highest average overlap score
        '''

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
        '''Builds a sequence of a path.

        :param path: Path
        :returns: Sequence built according to provided path
        '''

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
        '''Loads used sequences into graph nodes.

        :param path: FASTA/FASTQ file with sequences
        :param nodes: Graph nodes
        :param used_nodes: List of used nodes
        '''

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
        '''Determines file type.

        Method expects that provided file ends either with 'fasta'
        or 'fastq' extension.

        :param path: File
        :returns: File type according to file extension
        '''

        return 'fastq' if path.endswith('.fastq') else 'fasta'

    def prepare_test_data(self, path):
        '''Prepares test data according to final path.

        :param path: Path
        :returns: Corresponding TestData instance
        '''

        contigs = self.filter_contigs(path)
        longest_contig = max(contigs, key=lambda c: c.len)
        return TestData(len(contigs), longest_contig.len)

    @staticmethod
    def sequence_identity(ol: pafpy.PafRecord):
        '''Calculates sequence identity.

        :param ol: Overlap record
        :returns: Sequence identity of an overlap
        '''

        return ol.mlen / ol.blen

    @staticmethod
    def average_overlap_length(ol):
        '''Calculates average overlap length.

        :param ol: Overlap record
        :returns: Average overlap length
        '''

        return (ol.qend-ol.qstart+ol.tend-ol.tstart) / 2

    @staticmethod
    def overlap_score(ol):
        '''Calculates overlap score.

        :param ol: Overlap record
        :returns: Overlap score
        '''

        if isinstance(ol, pafpy.PafRecord):
            seq_id = Graph.sequence_identity(ol)
        else:
            seq_id = ol.seq_id
        return seq_id * Graph.average_overlap_length(ol)

    @staticmethod
    def extension_score(ol: Overlap, q, t):
        '''Calculates extension score.

        :param ol: Overlap record
        :param q: Extended sequence
        :param t: Extending sequence
        :returns: Extension score
        '''

        ol_score = Graph.overlap_score(ol)
        oh_1 = len(q) - ol.qend
        oh_2 = ol.tstart
        ext_len = len(t) - ol.tend
        return ol_score + ext_len/2 - (oh_1 + oh_2)/2

    @staticmethod
    def should_add_overlap(ol, read, compl_read, contigs, reads):
        '''Determines whether provided overlap should be added to graph.

        :param ol: Overlap record
        :param read: Read that is included in overlap
        :param compl_read: Complemented read
        :param contigs: Contigs
        :param reads: Reads
        :returns: True if overlap should be added to graph
        '''

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
