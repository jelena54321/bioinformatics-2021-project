class Overlap:

    def __init__(self, qstart, qend, tstart, tend, seq_id):
        self.qstart = qstart
        self.qend = qend

        self.tstart = tstart
        self.tend = tend

        self.seq_id = seq_id

class Node:

    COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, id, len):
        self.id = id
        self.len = len
        self.seq = None
        self.nodes = []
        self.overlaps = []

    def __str__(self):
        string = f'{self.id}: {len(self.nodes)}'
        for node in self.nodes:
            string += f'\n\t{node.id}'

        return string

    def __len__(self):
        return self.len

    def add_overlap(self, node, overlap):
        self.nodes.append(node)
        self.overlaps.append(overlap)

    def remove_overlap(self, node):
        idx = self.nodes.index(node)
        self.nodes.pop(idx)
        self.overlaps.pop(idx)

    def overlap_for_node(self, node):
        return self.overlaps[self.nodes.index(node)]

    @staticmethod
    def complement_id(id):
        idx = id.find('~')
        return '~' + id if idx == -1 else id[idx + 1:]

    @staticmethod
    def complement_sequence(seq):
        if seq is None: return None

        complemented_seq = []
        for base in seq:
            if base not in Node.COMPLEMENTS: complement = base
            else: complement = Node.COMPLEMENTS[base]

            complemented_seq.append(complement)

        return (''.join(complemented_seq))[::-1]

    @staticmethod
    def complement(node):
        return Node(Node.complement_id(node.id), len(node))
