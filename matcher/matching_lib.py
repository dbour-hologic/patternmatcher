""" MAIN LOGIC FOR STRING MATCHING
Author: David Bour
Date: 2016-06-27

### String Matching: Tuned for DNA ###
Naive implementation of string matching.

"""

def dna_is_a_match(pattern, subject_dna, start_index, mismatch_tolerance=0):

    """
    Function can be used for RNA or DNA; can also be used with IUPAC bases such as 
    'N','V','D' (etc). This is a naive implementation of string matching. 

    Args
        pattern - oligonucleotide sequence to match (str)
        subject_dna - subject sequence to match against pattern (str)
        start_index - where to begin searching (int)
        mismatch_tolerance - tolerance of mismatches, default is perfect match (int)
    Returns
        boolean. If match was found within tolerance level
        True - Found a match within the tolerance level
        False - Did not find a match within the tolerance level
    """

    # IUPAC Nucleotide Codes
    IUPAC = {
        "A" : ["A"], 
        "C" : ["C"], 
        "G" : ["G"],        
        "T" : ["T"],        
        "U" : ["U"],        
        "R" : ["G", "A"], 
        "Y" : ["T", "C"], 
        "K" : ["G", "T"], 
        "M" : ["A", "C"], 
        "S" : ["G", "C"], 
        "W" : ["A", "T"], 
        "B" : ["C", "G", "T"], 
        "D" : ["A", "G", "T"], 
        "H" : ["A", "C", "T"], 
        "V" : ["A", "C", "G"], 
        "N" : ["A", "C", "G", "T"]
    }

    #counter - a mismatch is defined as not being equivalent to any IUPAC code
    mismatch_found = 0

    if (start_index + len(pattern) > len(subject_dna)):
        return False

    for index in range(len(pattern)):
        
        if subject_dna[start_index+index] not in IUPAC.get(pattern[index]):

            mismatch_found += 1

            if mismatch_found > mismatch_tolerance:
                return False

    return True

def dna_first_match(pattern, subject_dna, mismatch_tolerance=0):
    return 10

def dna_match_all(pattern, subject_dna, mismatch_tolerance=0):
    return [0,1,2]



