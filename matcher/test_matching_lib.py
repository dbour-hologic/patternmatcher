""" 
Unit Test for matching_lib.py
Date: 2016-06-27
"""

import unittest
from matching_lib import dna_is_a_match, dna_first_match, dna_match_all 

class TestMatchMethods(unittest.TestCase):

	def test_dna_is_a_match(self):
		self.assertTrue(dna_is_a_match("ATCC","ATCCCC",0, 0), "Pattern should match with up to 0 mismatches to the subject.")
		self.assertTrue(dna_is_a_match("GTCC","ATCCCC",0, 1), "Pattern should match with up to 1 mismatches to the subject.")
		self.assertTrue(dna_is_a_match("ATGG","ATCCCC",0, 2), "Pattern should match with up to 2 mismatches to the subject.")
		self.assertTrue(dna_is_a_match("GGCC","ATCCCC",0, 3), "Pattern should match with up to 3 mismatches to the subject.")		
		self.assertTrue(dna_is_a_match("GGGG","ATCCCC",0, 4), "Pattern should match with up to 4 mismatches to the subject.")

		self.assertTrue(dna_is_a_match("NNNN","ATCCCC",0, 0), "Pattern should be able to take IUPAC codes.")
		self.assertTrue(dna_is_a_match("ATCN","ATCCCC",0, 0), "Pattern should be able to take IUPAC codes.")
		self.assertTrue(dna_is_a_match("BTCN","ATCCCC",0, 1), "Pattern should be able to take IUPAC codes.")

		self.assertTrue(dna_is_a_match("ATCC","ATCCCC",0,20), "Pattern should not fail even though threshold > length of pattern")

		self.assertFalse(dna_is_a_match("ATCC","ATCCCC",1,0), "Pattern should fail at the index shift.")

		self.assertFalse(dna_is_a_match("AAAAAAA","TA",0, 0), "Pattern should fail to match if its greater than the subject.")

	def test_dna_first_match(self):
		find_index = dna_first_match("ATCG","GATCGCCC", 0)
		self.assertEqual(find_index, 1, "First match should be at index 1 with 0 mismatches.")

		find_index_2 = dna_first_match("GTCG","GATCGCCC", 1)
		self.assertEqual(find_index_2, 1, "First match should at at index 1 with 1 mismatch.")

		find_index_3 = dna_first_match("ATCG","GGGGGATCG", 0)
		self.assertEqual(find_index_3, 5, "First match should be at index 5 with 0 mismatches.")

	def test_dna_match_all(self):
		list_of_indexes = dna_match_all("ATCG","ATCGATCGATCG", 0)
		self.assertEqual(list_of_indexes, [0,4,8], "Patterns should be found at position 0, 4, and 8 for 0 mismatches.")

		list_of_indexes_2 = dna_match_all("ATCG","ATCGGGGGATCG", 0)
		self.assertEqual(list_of_indexes_2, [0,8], "Patterns should be found at position 0 and 8 for 0 mismatches..")

		list_of_indexes_3 = dna_match_all("","ATCGGGGAGATCA", 0)
		self.assertEqual(list_of_indexes_3, [], "No matches should be found for blank patterns.")

		list_of_indexes_4 = dna_match_all("AT","ATATATAG", 1)
		self.assertEqual(list_of_indexes_4, [0, 2, 4, 6], "Patterns should be found at position 0, 2, 4, and 6 for 1 mismatch.")

if __name__ == '__main__':
	unittest.main()