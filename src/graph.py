from typing import List
import pafpy
from Bio import SeqIO, Seq, SeqRecord
from node import Node, Overlap
import random

SEQ_ID_MIN = 0.85
LEN_DELTA = 100

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
    def get_best_paths(self, l, num_elementa, next, current_path):
        num_elementa_in_next = 0

        for i in l:
            if i in current_path or i in next:
                continue
            next.append(i.node)
            num_elementa_in_next += 1
            if num_elementa_in_next == num_elementa:
                break

    def dfs(self, current_node: Node, current_path: List, all_found_paths: List, heuristic_depth: int):
        #dfs dosao do kraj, ili je zbog heuristike ili nea vise djece
        if heuristic_depth == 0 or len(current_node.nodes) == 0:
            #DEBUG
            print("grap.generate_paths.dfs >> FOUND ONE PATH :D")
            all_found_paths.append(current_path)
            return 

        # dobi sve elemente
        all_overlap_scores = []
        all_extension_scores = []

        for i in range(len(current_node.nodes)):
            all_overlap_scores.append(self.helper_sort(current_node.nodes[i], Graph.overlap_score(current_node.overlaps[i])))
            all_extension_scores.append(self.helper_sort(current_node.nodes[i], Graph.extension_score(current_node.overlaps[i])))

        sorted(all_overlap_scores, key=lambda el: el.score, reverse=True)
        sorted(all_extension_scores, key=lambda el: el.score, reverse=True)

        #Kreiramo next
        #koliko kojih elemenata ce uzimati za svaki od 3 nacina
        num_element = 20
        next = []

        self.get_best_paths(all_overlap_scores, num_element, next, current_path)
        self.get_best_paths(all_extension_scores, num_element, next, current_path)

        #monte_carlo
        num_elementa_in_next = 0
        len_of_all_elements = len(current_node.nodes)
        while(len_of_all_elements and num_elementa_in_next < num_element):
            len_of_all_elements -= 1
            el  = random.choice(current_node.nodes)
            if el in current_path or el in next:
                continue
            next.append(el)
            num_elementa_in_next += 1
       
        for i in next:
            current_path.append(i)
            self.dfs(i, current_path, all_found_paths, heuristic_depth - 1)
            #backtracking
            current_path.pop()
        

    def generate_paths(self):
        """
        TODO: Traverse through graph using depth first search and return list of found
        paths (list of node lists).
        """
        all_found_paths = []
        current_path = []
        heuristic_depth = 100 # trebamo Toonati :) 

        print("grap.generate_paths >> Finding_Paths")
        first_node = random.choice(list(self.contigs.values()))
        while(len(first_node.nodes) == 0):
             first_node = random.choice(list(self.contigs.values()))

        current_path.append(first_node)
        self.dfs(first_node, current_path, all_found_paths, heuristic_depth)
        
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

                if ol.strand == pafpy.Strand.Forward:
                    if t_left > q_left:
                        overlap = Overlap(
                            ol.tstart, ol.tend, ol.tlen,
                            ol.qstart, ol.qend, ol.qlen,
                            ol.mlen, ol.blen
                        )
                        t.add_overlap(q, overlap)

                        overlap = Overlap(
                            ol.qlen - ol.qend, ol.qlen - ol.qstart, ol.qlen,
                            ol.tlen - ol.tend, ol.tlen - ol.tstart, ol.tlen,
                            ol.mlen, ol.blen
                        )
                        compl_q.add_overlap(compl_t, overlap)

                    else:
                        overlap = Overlap(
                            ol.qstart, ol.qend, ol.qlen,
                            ol.tstart, ol.tend, ol.tlen,
                            ol.mlen, ol.blen
                        )
                        q.add_overlap(t, overlap)

                        overlap = Overlap(
                            ol.tlen - ol.tend, ol.tlen - ol.qstart, ol.tlen,
                            ol.qlen - ol.qend, ol.qlen - ol.qstart, ol.qlen,
                            ol.mlen, ol.blen
                        )
                        compl_t.add_overlap(compl_q, overlap)

                else:
                    if t_left > q_right:
                        overlap = Overlap(
                            ol.qstart, ol.qend, ol.qlen,
                            ol.tlen - ol.tend, ol.tlen - ol.tstart, ol.tlen,
                            ol.mlen, ol.blen
                        )
                        q.add_overlap(compl_t, overlap)

                        overlap = Overlap(
                            ol.tstart, ol.tend, ol.tlen,
                            ol.qlen - ol.qend, ol.qlen - ol.qstart, ol.qlen,
                            ol.mlen, ol.blen
                        )
                        t.add_overlap(compl_q, overlap)

                    else:
                        overlap = Overlap(
                            ol.tlen - ol.tend, ol.tlen - ol.tstart, ol.tlen,
                            ol.qstart, ol.qend, ol.qlen,
                            ol.mlen, ol.blen
                        )
                        compl_t.add_overlap(q, overlap)

                        overlap = Overlap(
                            ol.qlen - ol.qend, ol.qlen - ol.qstart, ol.qlen,
                            ol.tstart, ol.tend, ol.tlen,
                            ol.mlen, ol.blen
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
    def sequence_identity(ol):
        return ol.mlen / ol.blen

    @staticmethod
    def overlap_score(ol):
        avg_ol_len = (ol.qend - ol.qstart + ol.tend - ol.tstart) / 2
        return Graph.sequence_identity(ol) * avg_ol_len

    @staticmethod
    def extension_score(ol):
        ol_score = Graph.overlap_score(ol)
        oh_1 = ol.qlen - ol.qend
        oh_2 = ol.tstart
        ext_len = ol.tlen - ol.tend
        return ol_score + ext_len / 2 - (oh_1 + oh_2) / 2

    @staticmethod
    def should_add_overlap(ol, read, compl_read, contigs, reads):
        if len(read) == 0 and len(compl_read) == 0: return True

        r = read if len(read) > 0 else compl_read
        old_contig = r.nodes[0]
        old_overlap = r.overlaps[0]

        if Graph.overlap_score(old_overlap) > Graph.overlap_score(ol): return False

        r.remove_overlap(old_contig)

        compl_r = reads[Node.complement_id(r.id)]
        compl_c = contigs[Node.complement_id(old_contig.id)]
        compl_c.remove_overlap(compl_r)

        return True
