from collections import namedtuple


class SortHelper:
    '''Class used for node sorting.

    :param node: Node
    :param node: Node score
    '''

    def __init__(self, node, score):
        self.node = node
        self.score = score


class SearchState:
    '''Class that represents a search state.
    It is used for graph traversing.

    :param nodes: Expanded graph nodes
    :param path: Active node path
    :path len: Length of an active path
    :gap_len: Length of a filled gap between two contigs
    '''

    def __init__(self, nodes, path, path_len, gap_len):
        self.nodes = nodes
        self.path = path
        self.path_len = path_len
        self.gap_len = gap_len


'''Tuple used for testing.

:param n_contigs: Number of contigs in the path
:param contig_len: Length of the longest included contig
'''
TestData = namedtuple('TestData', ['n_contigs', 'contig_len'])
