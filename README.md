# Improvement of partially assembled genome using long reads
Program based on the method [Highly Efficient Repeat Assembly (HERA)](https://www.biorxiv.org/content/10.1101/345983v1) that uses long reads to connect contigs to longer sequences.

## Dependencies
- python3.8-dev
- python-virtualenv

## Installation
```bash
git clone https://github.com/jelena54321/bioinformatics-2021-project.git
cd bioinformatics-2021-project
make install
```

## Usage
Activate virtual environment:
```bash
cd bioinformatics-2021-project
. ./venv/bin/activate
```

Run program:
```
    python3 ./src/main.py <contigs> <reads> <reads_to_contigs> <reads_to_reads> <out>

        <contigs>
            Contigs in FASTA/FASTQ format
        <reads>
            Reads in FASTA/FASTQ format
        <reads_to_contigs>
            <reads> aligned to <contigs> in PAF format
        <reads_to_reads>
            <reads> aligned to <reads> in PAF format
        <out>
            Output file in FASTA format.
```

## Testing

## Clean
Clean all data and delete virtual environment:
```bash
make clean
```
