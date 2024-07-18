#!/usr/bin/env python

import math
import pandas as pd
import sys
from scipy.signal import find_peaks
import json
import numpy as np
import pybedtools




"""
Read bedfile into pandas dataframe
"""
def readBed(file):
	df = pd.read_csv(file, sep='\t', header=None)
	header = ['chrom', 'chromStart', 'chromEnd', 'score']
	df.columns = header[:len(df.columns)]
	return df


## based on termseq_peaks
class BedObject:
	def __init__(self, bedfile, strand="+"):
		self.df = readBed(bedfile)
		self.peaks_by_chrom = {}
		self.meta_by_chrom = {}
		self.starts_by_chrom = {}
		self.chroms = []
		self.strand = strand

	def call_peaks(self):

		for chrom in self.df['chrom'].unique():
			idx = self.df['chrom'] == chrom
			starts = self.df.loc[idx, "chromStart"]
			x = self.df.loc[idx, "score"]
			## optional: expand to include each position
			peaks, meta = find_peaks(x, prominence=(None, None), width=(1, None), rel_height=0.75)
			# meta keeps the prominences, plateau_size etc.
			self.peaks_by_chrom[chrom] = peaks
			self.meta_by_chrom[chrom] = meta
			self.starts_by_chrom[chrom] = starts
			self.chroms.append(chrom)

		return self


	# from termseq_peaks
	def peaks_to_bed(self):
		"""
		Call peaks into the internal format, and then output as a narrowPeak file.

		Returns
		-------
		pybedtools.BedTool object sorted by score
		"""

		def gen():
			for chrom in self.chroms:
				starts = self.starts_by_chrom[chrom]
				left_ips = self.meta_by_chrom[chrom]["left_ips"]
				right_ips = self.meta_by_chrom[chrom]["right_ips"]
				left_bases = self.meta_by_chrom[chrom]["left_bases"]
				right_bases = self.meta_by_chrom[chrom]["right_bases"]
				prominences = self.meta_by_chrom[chrom]["prominences"]
				widths = self.meta_by_chrom[chrom]["widths"]
				peaks = self.peaks_by_chrom[chrom]

				xp = np.arange(len(starts))

				ileft_ips = np.interp(left_ips, xp, starts).round().astype(int)
				iright_ips = np.interp(right_ips, xp, starts).round().astype(int)
				ipeaks = np.interp(peaks, xp, starts).round().astype(int)

				# positions where left < right (i.e., not from end of chrom to beginning)
				idx = ileft_ips <= iright_ips
				ileft_ips = ileft_ips[idx]
				iright_ips = iright_ips[idx]
				ipeaks = ipeaks[idx]
				widths = widths[idx]
				prominences = prominences[idx]

				# number of removed peaks is the sum of False values
				n_removed = sum(~idx)
				if n_removed:
                			print("Peaks removed due to start/stop problems: {0}".format(n_removed))


				for start, stop, peak, prominence, width in zip(ileft_ips, iright_ips, ipeaks, prominences, widths):

					# This uses the promience as the score.
					p = str(prominence)

					yield pybedtools.create_interval_from_list(
						[chrom, str(start), str(stop), ".", p, self.strand, p, "-1", "-1", str(peak - start)])

		# Ensure we're coord-sorted for the merging step
		x = pybedtools.BedTool(gen()).sort()
		x = merge_narrowbed(x, self.strand)

		# But the output needs to be sorted by score
		return sort_by_score(x)


# from termseq-peaks
def merge_narrowbed(peaks, strand, additional_kwargs={"d": 1, "o": "max"}):
	x = (peaks.cut([0, 1, 2, 4]).merge(c=4, **additional_kwargs)).to_dataframe()
	x["score"] = "."
	x["strand"] = strand
	y = pybedtools.BedTool.from_dataframe(x)
	return y.each(to_narrowpeak).saveas()

# from termseq-peaks
def to_narrowpeak(f):
	"""
	Convert a feature into narrowPeak format, with signal and pval equivalent
	to the score.
	"""
	return pybedtools.create_interval_from_list([ f.chrom, str(f.start), str(f.stop), ".", f.name, f.strand, f.name, f.name, "-1",
		str(int((f.stop - f.start) / 2))])

# from termseq-peaks
def sort_by_score(x):
	"""
	Sort a BedTool object by the score column.
	"""
	df = pybedtools.BedTool(x).to_dataframe()
	df = df.sort_values("score", ascending=False)
	return pybedtools.BedTool.from_dataframe(df)


if __name__ == "__main__":

	bedfile = sys.argv[1]
	outputfile = sys.argv[2]
	strand = sys.argv[3]
	
	bed = BedObject(bedfile, strand)
	bed.call_peaks().peaks_to_bed().moveto(outputfile)
