import unittest
import dnaSequence

class TestDnaSequence(unittest.TestCase):
    def test_neighbors(self):
        sequence = "ACG"
        expected = ["CCG", "TCG", "GCG", "AAG", "ATG", "AGG", "ACA", "ACC", "ACT", "ACG"]
        actual = dnaSequence.neighbors(sequence, 1)
        self.assertEqual(set(expected), set(actual))

    def test_contains(self):
        sequence = "ATTTTGGC"
        tofind = "ATTT"
        self.assertEqual(True, dnaSequence.contains(sequence, tofind, 0))

        tofind = "AGTT"
        self.assertEqual(True, dnaSequence.contains(sequence, tofind, 1))

        tofind = "AGG"
        self.assertEqual(False, dnaSequence.contains(sequence, tofind, 0))

        tofind = "AGGT"
        self.assertEqual(False, dnaSequence.contains(sequence, tofind, 1))

        tofind = "C"
        self.assertEqual(True, dnaSequence.contains(sequence, tofind, 0))

    def test_isMotif(self):
        sequences = ["GATG", "GATG"]
        tofind = "AT"
        self.assertEqual(True, dnaSequence.isMotif(sequences, tofind, 0))

        sequences = ["GAT", "GAG"]
        toFind = "AT"
        self.assertEqual(True, dnaSequence.isMotif(sequences, toFind, 1))

        sequences = ["GAT", "GAG"]
        toFind = "AT"
        self.assertEqual(False, dnaSequence.isMotif(sequences, toFind, 0))

    def test_findMotifs(self):
        sequences = ["AT", "GA"]
        expected = set(["AA", "GT"])
        actual = dnaSequence.findMotifs(sequences, 2, 1)
        self.assertEqual(expected, actual)

        sequences = ["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"]
        expected = set(["ATA", "ATT", "GTT", "TTT"])
        actual = dnaSequence.findMotifs(sequences, 3, 1)
        self.assertEqual(expected, actual)

        sequences = ["AAAAA", "AAAAA","AAAAA"]
        expected = set(["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"])
        actual = dnaSequence.findMotifs(sequences, 3, 3)
        self.assertEqual(expected, actual)

        sequences = ["ACGT", "ACGT", "ACGT"]
        expected = set(["ACG", "CGT"])
        actual = dnaSequence.findMotifs(sequences, 3, 0)
        self.assertEqual(expected, actual)

        sequences = ["AAAAA", "AAAAA", "AACAA"]
        expected = set([])
        actual = dnaSequence.findMotifs(sequences, 3, 0)
        self.assertEqual(expected, actual)

        sequences = ["AACAA", "AAAAA", "AAAAA"]
        expected = set([])
        actual = dnaSequence.findMotifs(sequences, 3, 0)
        self.assertEqual(expected, actual)

        sequences = ["AAAAA", "AAAAA", "AAAAA"]
        expected = set(["AAA", "AAC", "AAG", "AAT", "ACA", "AGA", "ATA", "CAA", "GAA", "TAA"])
        actual = dnaSequence.findMotifs(sequences, 3, 1)
        self.assertEqual(expected, actual)

    def test_subs(self):
        string = ["GATATATGCATATACTT"]
        sub = ["ATAT"]

        expected = [2,4,10]
        actual = dnaSequence.subs(string, sub)
        self.assertEqual(expected, actual)
