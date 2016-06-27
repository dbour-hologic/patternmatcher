""" 
Unit Test for matching_lib.py
Date: 2016-06-27
"""

import unittest
from matching_lib import dna_is_a_match, dna_first_match, dna_match_all 

class TestMatchMethods(unittest.TestCase):

	def test_dna_is_a_match(self):
		self.assertTrue(dna_is_a_match("ATCC","ATCCCC",0, 0), "Pattern should match with 0 mismatches to the subject.")
		self.assertTrue(dna_is_a_match("GTCC","ATCCCC",0, 1), "Pattern should match with 1 mismatches to the subject.")
		self.assertTrue(dna_is_a_match("ATGG","ATCCCC",0, 2), "Pattern should match with 2 mismatches to the subject.")
		self.assertTrue(dna_is_a_match("GGCC","ATCCCC",0, 3), "Pattern should match with 3 mismatches to the subject.")		
		self.assertTrue(dna_is_a_match("GGGG","ATCCCC",0, 4), "Pattern should match with 4 mismatches to the subject.")
		self.assertTrue(dna_is_a_match("ATCC","ATCCCC",0,20), "Pattern should not fail even though threshold > length of pattern")
		self.assertFalse(dna_is_a_match("ATCC","ATCCCC",1,0), "Pattern should fail at the index shift.")


	def test_dna_first_match(self):
		find_index = dna_first_match("ATCG","GATCGCCC")
		self.assertEqual(find_index, 1, 0)

	def test_dna_match_all(self):
		list_of_index = dna_match_all("ATCG","ATCGATCGATCG")
		self.assertEqual(list_of_index, [0,4,8])

if __name__ == '__main__':
	unittest.main()