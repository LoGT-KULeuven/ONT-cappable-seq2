#!/usr/bin/env python3

from Bio import SeqIO
import sys
import os


def removeEmptySequences(infile, outfile):
        record_iterator = SeqIO.parse(infile, "fastq")

        with open(outfile, "w") as out_ok:
            for record in record_iterator:
                if len(str(record.seq)) != 0:
                    SeqIO.write(record, out_ok, 'fastq')


if __name__ == "__main__":
	infile = sys.argv[1]
	outfile = sys.argv[2]
	removeEmptySequences(infile, outfile)
