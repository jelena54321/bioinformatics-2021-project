from collections import namedtuple


class SortHelper:

    def __init__(self, node, score):
        self.node = node
        self.score = score


class SearchState:

    def __init__(self, nodes, path, path_len, gap_len):
        self.nodes = nodes
        self.path = path
        self.path_len = path_len
        self.gap_len = gap_len


TestData = namedtuple('TestData', ['n_contigs', 'contig_len'])
