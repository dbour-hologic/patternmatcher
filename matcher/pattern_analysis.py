""" Pattern Analysis Program
Author: David Bour
Date: 2016-06-28

"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from matching_lib import dna_match_all

class RunIdentifier(object):

	""" 
	The RunIdentifier class is used to group the results of 
	hits to a specific pattern; a one to one relationship. 
	"""

	def __init__(self, id, subject, pattern, pattern_id, hit_locations, total_hits, mismatch_tolerance_set):

		""" Constructor for RunIdentifier

		Args:
			id - the FASTA parsed ID of the subject (str)
			subject - the sequence to be queried (str)
			pattern - the pattern to be queried against the subject (str)
			pattern_id - the FASTA parsed ID of the pattern (str)
			hit_locations - all starting positions of matches (list)
			total_hits - total number of times the specific pattern hit per sequence (int)
			mismatch_tolerance_set - the mismatch tolerance set for this particular run

		"""

		self.id = id
		self.subject = subject
		self.pattern = pattern
		self.pattern_id = pattern_id
		self.hit_locations = hit_locations
		self.total_hits = total_hits
		self.mismatch_tolerance_set = mismatch_tolerance_set

	def print_friendly(self):

		for location in self.hit_locations:
			end_point = len(self.pattern) + location - 1
			print "Pattern: {0}".format(self.pattern)
			print "Subject: {0}".format(self.subject[location:end_point+1])
			print "Start Location: {0} End Location: {1}".format(location+1, end_point)

	def print_more(self):

		for location in self.hit_locations:

			end_point = len(self.pattern) + location - 1

			extend_right, extend_left = 4, 4

			if location-4 <= 0:
				extend_left = location
			if end_point+4 > len(self.subject):
				extend_right = len(self.subject) - end_point


			pattern_format = " " * extend_left + " " + self.pattern + " " + " " * extend_right
			subject_format = self.subject[0:extend_left] + " " + self.subject[location:end_point+1] + " " + self.subject[end_point:end_point+extend_right]

			print "Pattern: {0}".format(pattern_format)
			print "Subject: {0}".format(subject_format)


class PatternAnalysis(object):

	""" Main Analysis Component

	Makes use of Biopython as the FASTA parser and a custom DNA matching
	library for pattern matching.

	"""

	def __init__(self):
		"""
		Constructor for Pattern Analysis Program

		list_of_queries (dict) - Holds a unique sequence ID (key) along with a list
		of RunIdentifier objects (value). Helps to keep all of the queries against
		the subjects.

		loaded_patterns (list) - Holds all of the patterns in a SeqRecord object format

		loaded_subjects (list) - Holds all of the target sequences in a SeqRecord object format

		"""
		self.list_of_queries = {}
		self.loaded_patterns = []
		self.loaded_subjects = []

	def _load_patterns(self, pattern_file):
		"""
		Reads a .txt file to load patterns to query against
		the subjects. Uses Biopython module for FASTA parsing.

		Args:
			pattern_file - directory path to pattern file (str)

		"""
		try:
			for pat_records in SeqIO.parse(pattern_file, "fasta"):
				self.loaded_patterns.append(pat_records)
		except IOError:
			print "Failed to find pattern file at '{0}'".format(pattern_file)

	def _load_subjects(self, subject_file):
		"""
		Reads a .txt file to load subject sequences to query against
		the patterns. Uses Biopython module for FASTA parsing.

		Args:
			subject_file - directory path to pattern file (str)

		"""
		try:
			for records in SeqIO.parse(subject_file, "fasta"):
				self.loaded_subjects.append(records)
		except IOError:
			print "Failed to find subject file at '{0}'".format(subject_file)

	def run(self, pattern_file, subject_file, mismatch_tolerance=0):
		"""
		CORE LOGIC OF PROGRAM
		> Sets up the environment and runs the program.

		Args
			pattern_file - file path of the pattern file in FASTA <*.fasta> format (str)
			subject_file - file path of the subject file in FASTA <*.fasta> format (str)
			mismatch_tolerance - the amount of mismatches to tolerate (int)
		"""

		# Load up the data from the files
		self._load_patterns(pattern_file)
		self._load_subjects(subject_file)

		# Run each subject against a set of oligonucleotide patterns
		for subject in self.loaded_subjects:

			# The sequence ID
			run_id = subject.id
			# The sequence itself
			run_seq = str(subject.seq).rstrip()
			# Create a placeholder to hold all patterns searched against this
			# particular sequence.
			self.list_of_queries[run_id] = []

			for pattern in self.loaded_patterns:
				# The pattern ID
				pattern_id = pattern.id
				# The pattern sequence
				pattern_seq = str(pattern.seq).rstrip()
				# Core matching algorithm
				results = dna_match_all(pattern_seq, run_seq, mismatch_tolerance)
				# Create the sequence-pattern pair object to store
				run_create = RunIdentifier(
											run_id, 
											run_seq, 
											pattern_seq, 
											pattern_id, 
											results, 
											len(results), 
											mismatch_tolerance
										  )
										
				self.list_of_queries[run_id].append(run_create)

	def size(self):
		return len(self.list_of_queries)

	def print_features(self):
		for key, value in self.list_of_queries.iteritems():
			print value[0].print_more()
			print "________"

if __name__ == '__main__':

	test_run = PatternAnalysis()
	test_run.run('data/patterns.txt', 'data/2016FASTA.fasta', 3)
	test_run.print_features()




