class Overlap:

    def __init__(self, qstart, qend, tstart, tend, seq_id):
        self.qstart = qstart
        self.qend = qend
        self.tstart = tstart
        self.tend = tend
        self.seq_id = seq_id


class Node:

    def __init__(self, id, len):
        self.id = id
        self.len = len
        self.seq = None
        self.nodes = []
        self.overlaps = []

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
        return '~' + id if idx == -1 else id[idx+1:]

    @staticmethod
    def complement_sequence(seq):
        return None if seq is None else seq.reverse_complement()

    @staticmethod
    def complement(node):
        return Node(Node.complement_id(node.id), len(node))
