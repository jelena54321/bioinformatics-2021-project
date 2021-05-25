import argparse
import sys
sys.path.append('../src')
from graph import Graph


def prepare_test_data_string(test_data):
    return (
        '\nTest Data:\n'
        f'* number of contigs in the result path: {test_data.n_contigs}\n'
        f'* length of the longest included contig: {test_data.contig_len}\n'
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('contigs', type=str)
    parser.add_argument('reads', type=str)
    parser.add_argument('reads_to_contigs', type=str)
    parser.add_argument('reads_to_reads', type=str)
    parser.add_argument('out', type=str)
    args = parser.parse_args()

    graph = Graph.construct(args.reads_to_contigs, args.reads_to_reads)
    paths, path_lens = graph.generate_paths()
    test_data = graph.generate_sequence(
        paths, path_lens,
        args.contigs, args.reads, args.out
    )

    print(prepare_test_data_string(test_data))


if __name__ == '__main__':
    main()
