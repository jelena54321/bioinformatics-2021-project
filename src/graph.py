import pafpy
from Bio import SeqIO
from Bio import pairwise2

SEQ_ID_MIN = 0.85

class Node:

    COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, id, seq=None):
        self.id = id
        self.seq = seq
        self.overlaps = []
        self.ol_scores = []
        self.ext_scores = []

    def __str__(self):
        string = f'{self.id}: {len(self)}'
        for overlap in self.overlaps:
            string += f'\n\t+{overlap.id}'

        return string

    def __len__(self):
        return len(self.overlaps)

    def add_overlap(self, node, ol_score, ext_score):
        self.overlaps.append(node)
        self.ol_scores.append(ol_score)
        self.ext_scores.append(ext_score)

    def remove_overlap(self, node):
        if node not in self.overlaps: return

        idx = self.overlaps.index(node)
        self.overlaps.pop(idx)
        self.ol_scores.pop(idx)
        self.ext_scores.pop(idx)

    def contains_overlap(self, node):
        return node in self.overlaps

    def complement_id(self):
        idx = self.id.find('~')
        return '~' + self.id if idx == -1 else self.id[idx + 1:]

    @staticmethod
    def complement(id, seq=None):
        # TODO: CHECK
        # if it is better to save complemented sequence rather than
        # calculate every time: time vs. memory
        return Node(f'~{id}', Node.__complement_sequence(seq))

    @staticmethod
    def __complement_sequence(seq):
        if seq is None: return None

        complemented_seq = []
        for base in seq:
            if base not in Node.COMPLEMENTS: complement = base
            else: complement = Node.COMPLEMENTS[base]

            complemented_seq.append(complement)

        return complemented_seq

class Graph:

    def __init__(self, contigs, reads):
        self.contigs = contigs
        self.reads = reads

    @staticmethod
    def construct(contigs_path, reads_path, reads_to_contigs_path, reads_to_reads_path):
        contigs = dict()
        reads = dict()

        # TODO: CHECK
        # if we need to save sequences - takes a lot of time and memory
        """
        with open(contigs_path) as handle:
            for contig in SeqIO.parse(handle, 'fasta'):
                contigs[contig.id] = Node(contig.id, contig.seq)
                contigs[contig.id] = Node.complement(contig.id, contig.seq)
        print(f'graph.Graph.construct >> Loaded contigs.')

        with open(reads_path) as handle:
            for read in SeqIO.parse(handle, 'fasta'):
                reads[read.id] = Node(read.id, read.seq)
                reads[read.id] = Node.complement(read.id, read.seq)
        print(f'graph.Graph.construct >> Loaded reads.')
        """

        with pafpy.PafFile(reads_to_contigs_path) as paf:
            for align in paf:
                if align.is_unmapped(): continue
                if Graph.__sequence_identity(align) < SEQ_ID_MIN: continue

                if align.strand != pafpy.Strand.Reverse: continue
                if align.blen != align.mlen: continue

                if align.qname not in contigs:
                    contig = Node(align.qname)
                    contigs[align.qname] = contig
                    compl_contig = Node.complement(align.qname)
                    contigs[compl_contig.id] = compl_contig
                else:
                    contig = contigs[align.qname]
                    compl_contig = contigs[contig.complement_id()]

                if align.tname not in reads:
                    read = Node(align.tname)
                    reads[align.tname] = read
                    compl_read = Node.complement(align.tname)
                    reads[compl_read.id] = compl_read
                else:
                    read = reads[align.tname]
                    compl_read = reads[read.complement_id()]

                contig_left = align.qstart
                contig_right = align.qlen - align.qend
                read_left = align.tstart
                read_right = align.tlen - align.tend

                contig_contains_read = read_left <= contig_left and read_right <= contig_right
                if contig_contains_read: continue

                ol_score = Graph.__overlap_score(align)
                if len(read) > 0 or len(compl_read) > 0:
                    r = read if len(read) > 0 else compl_read
                    old_contig = r.overlaps[0]
                    old_ol_score = r.ol_scores[0]

                    if old_ol_score > ol_score:
                        continue
                    else:
                        r.remove_overlap(old_contig)

                        compl_r = reads[r.complement_id()]
                        compl_c = contigs[old_contig.complement_id()]
                        compl_c.remove_overlap(compl_r)

                # TODO: CHECK
                # start and end indices on target when reversed strand
                if align.strand == pafpy.Strand.Forward:
                    if read_left > contig_left:
                        ext_score = Graph.__extension_score(ol_score, contig_left, read_right, contig_right)
                        read.add_overlap(contig, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, contig_left, read_right, read_left)
                        compl_contig.add_overlap(compl_read, ol_score, ext_score)

                    else:
                        ext_score = Graph.__extension_score(ol_score, read_left, contig_right, read_right)
                        contig.add_overlap(read, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, read_left, contig_right, contig_left)
                        compl_read.add_overlap(compl_contig, ol_score, ext_score)

                else:
                    if read_left > contig_left:
                        ext_score = Graph.__extension_score(ol_score, contig_left, read_right, contig_right)
                        compl_read.add_overlap(contig, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, contig_left, read_right, read_left)
                        compl_contig.add_overlap(read, ol_score, ext_score)

                    else:
                        ext_score = Graph.__extension_score(ol_score, read_left, contig_right, read_right)
                        contig.add_overlap(compl_read, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, read_left, contig_right, contig_left)
                        read.add_overlap(compl_contig, ol_score, ext_score)

        print(f'graph.Graph.construct >> Generated overlaps between contigs and reads.')

        with pafpy.PafFile(reads_to_reads_path) as paf:
            for align in paf:
                if align.is_unmapped(): continue
                if Graph.__sequence_identity(align) < SEQ_ID_MIN: continue

                if align.qname not in reads:
                    q_read = Node(align.qname)
                    reads[align.qname] = q_read
                    compl_q_read = Node.complement(align.qname)
                    reads[compl_q_read.id] = compl_q_read
                else:
                    q_read = reads[align.qname]
                    compl_q_read = reads[q_read.complement_id()]

                if align.tname not in reads:
                    t_read = Node(align.tname)
                    reads[align.tname] = t_read
                    compl_t_read = Node.complement(align.tname)
                    reads[compl_t_read.id] = compl_t_read
                else:
                    t_read = reads[align.tname]
                    compl_t_read = reads[t_read.complement_id()]

                q_read_left = align.qstart
                q_read_right = align.qlen - align.qend
                t_read_left = align.tstart
                t_read_right = align.tlen - align.tend

                q_contains_t = q_read_left >= t_read_left and q_read_right >= q_read_right
                t_contains_q = q_read_left <= t_read_left and q_read_right <= t_read_right
                if q_contains_t or t_contains_q: continue

                ol_score = Graph.__overlap_score(align)
                if align.strand == pafpy.Strand.Forward:
                    if q_read_left > t_read_left:
                        if q_read.contains_overlap(t_read): continue

                        ext_score = Graph.__extension_score(ol_score, t_read_left, q_read_right, t_read_right)
                        q_read.add_overlap(t_read, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, t_read_left, q_read_right, q_read_left)
                        compl_t_read.add_overlap(compl_q_read, ol_score, ext_score)

                    else:
                        if t_read.contains_overlap(q_read): continue

                        ext_score = Graph.__extension_score(ol_score, q_read_left, t_read_right, q_read_right)
                        t_read.add_overlap(q_read, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, q_read_left, t_read_right, t_read_left)
                        compl_q_read.add_overlap(compl_t_read, ol_score, ext_score)

                else:
                    if q_read_left > t_read_left:
                        if q_read.contains_overlap(compl_t_read): continue

                        ext_score = Graph.__extension_score(ol_score, t_read_left, q_read_right, t_read_right)
                        q_read.add_overlap(compl_t_read, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, t_read_left, q_read_right, q_read_left)
                        t_read.add_overlap(compl_q_read, ol_score, ext_score)

                    else:
                        if compl_t_read.contains_overlap(q_read): continue

                        ext_score = Graph.__extension_score(ol_score, q_read_left, t_read_right, q_read_right)
                        compl_t_read.add_overlap(q_read, ol_score, ext_score)

                        ext_score = Graph.__extension_score(ol_score, q_read_left, t_read_right, t_read_left)
                        compl_q_read.add_overlap(t_read, ol_score, ext_score)

        print('graph.Graph.construct >> Generated overlaps between reads.')

        return Graph(contigs, reads)

    def generate_paths(self):
        pass

    def generate_sequences(self):
        pass

    @staticmethod
    def __sequence_identity(align):
        # TODO: CHECK
        # BLAST-identity: return align.mlen / align.blen
        # vs.
        # dv: approx. per-base divergence
        return 1 - align.get_tag('dv').value

    @staticmethod
    def __overlap_score(align):
        avg_ol_len = (align.query_aligned_length + align.target_aligned_length) / 2
        return Graph.__sequence_identity(align) * avg_ol_len

    @staticmethod
    def __extension_score(ol_score, oh_1, oh_2, ext_len):
        return ol_score + ext_len / 2 - (oh_1 + oh_2) / 2
