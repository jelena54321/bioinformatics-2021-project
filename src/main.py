import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('contigs', type=str)
    parser.add_argument('reads', type=str)
    parser.add_argument('reads_to_contigs', type=str)
    parser.add_argument('reads_to_reads', type=str)
    parser.add_argument('out', type=str)
    args = parser.parse_args()

if __name__ == '__main__':
    main()
