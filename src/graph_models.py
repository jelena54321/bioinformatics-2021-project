class Overlap:
    '''Class that represents a single overlap between two sequences.
    This class assumes that `t` always extends `q` and that both
    sequences are on the same strand.

    :param qstart: Overlap start of the extended sequence
    :param qend: Overlap end of the extended sequence
    :param tstart: Overlap start of the extending sequence
    :param tend: Overlap end of the extending sequence
    :param seq_id: Sequence identity of the overlap
    '''

    def __init__(self, qstart, qend, tstart, tend, seq_id):
        self.qstart = qstart
        self.qend = qend
        self.tstart = tstart
        self.tend = tend
        self.seq_id = seq_id


class Node:
    '''Class that represents a graph node, i.e. a sequence.

    :param id: Sequence ID
    :param len: Sequence length (in nucleobases)
    :param seq: Object that encapsulates sequence string
    :param nodes: List of nodes, i.e. extending sequences
    :param overlaps: List of overlaps
    '''

    def __init__(self, id, len):
        self.id = id
        self.len = len
        self.seq = None
        self.nodes = []
        self.overlaps = []

    def __len__(self):
        '''Gets sequence length.

        :returns: Sequence length
        '''

        return self.len

    def add_overlap(self, node, overlap):
        '''Adds new node and correponding overlap data.

        :param node: Node
        :param overlap: Overlap
        '''

        self.nodes.append(node)
        self.overlaps.append(overlap)

    def remove_overlap(self, node):
        '''Removes provided node and corresponding overlap data.

        :param node: Node
        '''

        idx = self.nodes.index(node)
        self.nodes.pop(idx)
        self.overlaps.pop(idx)

    def overlap_for_node(self, node):
        '''Gets overlap data for provided node.

        :param node: Node
        :returns: Corresponding overlap
        '''

        return self.overlaps[self.nodes.index(node)]

    @staticmethod
    def complement_id(id):
        '''Gets ID of a complemented node.

        :param id: Node ID
        :returns: ID of a complemented node
        '''

        idx = id.find('~')
        return '~' + id if idx == -1 else id[idx+1:]

    @staticmethod
    def complement_sequence(seq):
        '''Gets reverse complemented sequence.

        :param seq: Sequence
        :returns: Reverse complemented sequence
        '''

        return None if seq is None else seq.reverse_complement()

    @staticmethod
    def complement(node):
        '''Creates node with reverse complemented sequence.

        :param node: Node
        :returns: New reverse complemented node
        '''

        return Node(Node.complement_id(node.id), len(node))
