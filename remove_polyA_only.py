#!/usr/bin/env python3

from Bio import SeqIO
import sys
import os

def removepolyAonlyKeepBoth(file):
	record_iterator = SeqIO.parse(file, "fastq")

	filename = os.path.basename(file)
	folder = os.path.dirname(file)

	filename_ok = "passed_" + filename
	filename_notok = "removed_" +  filename

	with open(filename_ok, "w") as out_ok:
		with open(filename_notok, "w") as out_notok:

			for record in record_iterator:
				if len(str(record.seq).replace("A","")) != 0:
					SeqIO.write(record, out_ok, 'fastq')
				else:
					SeqIO.write(record, out_notok, 'fastq')

def removepolyAonly(infile, outfile):
        record_iterator = SeqIO.parse(infile, "fastq")

        with open(outfile, "w") as out_ok:
            for record in record_iterator:
                if len(str(record.seq)) != 0:
                    SeqIO.write(record, out_ok, 'fastq')


if __name__ == "__main__":
	infile = sys.argv[1]
	outfile = sys.argv[2]
	removepolyAonly(infile, outfile)
