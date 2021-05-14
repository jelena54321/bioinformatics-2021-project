import argparse
from graph import Graph


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('contigs', type=str)
    parser.add_argument('reads', type=str)
    parser.add_argument('reads_to_contigs', type=str)
    parser.add_argument('reads_to_reads', type=str)
    parser.add_argument('out', type=str)
    args = parser.parse_args()

    graph = Graph.construct(args.reads_to_contigs, args.reads_to_reads)
    paths = graph.generate_paths()
    graph.generate_sequence(paths, args.contigs, args.reads, args.out)


if __name__ == '__main__':
    main()
