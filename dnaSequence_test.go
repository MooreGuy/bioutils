package bioutils

import (
	"reflect"
	"testing"
)

/*
func TestCreateGenomeFromText(t *testing.T) {
	var actual, err = CreateDNASequence("ACGT")
	if err != nil {
		t.Error("Didn't expect an error from creating dna sequence.")
	}

	var expected = "ACGT"
	if expected != string(actual) {
		t.Error("Expected: ", expected, " does not match actual: ", actual)
	}

	actual, err = CreateDNASequence("ACGF")
	if err == nil {
		t.Error("Expected an error from giving an invalid base.")
	}
}
*/

func TestLeftPad(t *testing.T) {
	padding := dnaSequence([]nucleotide{nucleotide(PadNucleotide()), nucleotide(PadNucleotide()), nucleotide(PadNucleotide())})
	base, _ := CreateDNASequence("GGG")
	expected := append(padding, base...)
	actual := LeftPad(base, 3)
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

func TestHammingDistance(t *testing.T) {
	expectedDistance := 3
	sequence1, err := CreateDNASequence("GGGCCGTTGGT")
	if err != nil {
		t.Error("Didn't expect an error when creating a sequence")
	}

	sequence2, err := CreateDNASequence("GGACCGTTGAC")
	if err != nil {
		t.Error("Didn't expect an error when creating a sequence")
	}

	actualDistance := sequence1.HammingDistance(sequence2)
	if expectedDistance != actualDistance {
		t.Error("Expected: ", expectedDistance, "Actual: ", actualDistance)
	}
}

func TestGenerateReverseComplements(t *testing.T) {
	sequence, _ := CreateDNASequence("AT")
	expected, _ := CreateDNASequence("AT")
	actual := sequence.GenerateReverseComplement()
	if !reflect.DeepEqual(actual, expected) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}

	sequence, _ = CreateDNASequence("GGG")
	expected, _ = CreateDNASequence("CCC")
	actual = sequence.GenerateReverseComplement()
	if !reflect.DeepEqual(actual, expected) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}

	sequence, _ = CreateDNASequence("TGTC")
	expected, _ = CreateDNASequence("GACA")
	actual = sequence.GenerateReverseComplement()
	if !reflect.DeepEqual(actual, expected) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}

	sequence, _ = CreateDNASequence("TAT")
	expected, _ = CreateDNASequence("ATA")
	actual = sequence.GenerateReverseComplement()
	if !reflect.DeepEqual(actual, expected) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}
}

func TestFrequentKmersWithMismatchesAndReverseComplements(t *testing.T) {
	sequence, _ := CreateDNASequence("ACGTTGCATGTCGCATGATGCATGAGAGCT")
	expected, _ := DNASequences([]string{"ACAT", "ATGT"})
	actual := FrequentKmersWithMismatchesAndReverseComplements(sequence, 4, 1)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}

	sequence, _ = CreateDNASequence("AAT")
	expected, _ = DNASequences([]string{"AAT", "ATT"})
	actual = FrequentKmersWithMismatchesAndReverseComplements(sequence, 3, 0)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}
}

func TestFrequentKmersWithMismatches(t *testing.T) {
	sequence, _ := CreateDNASequence("ACGTTGCATGTCGCATGATGCATGAGAGCT")
	expected, _ := DNASequences([]string{"GATG", "ATGC", "ATGT"})
	actual := FrequentKmersWithMismatches(sequence, 4, 1)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}

	sequence, _ = CreateDNASequence("AAAAAAAAAA")
	expected, _ = DNASequences([]string{"AA", "AC", "AG", "CA", "AT", "GA", "TA"})
	actual = FrequentKmersWithMismatches(sequence, 2, 1)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}

	sequence, _ = CreateDNASequence("AGTCAGTC")
	expected, _ = DNASequences([]string{"TCTC", "CGGC", "AAGC", "TGTG", "GGCC",
		"AGGT", "ATCC", "ACTG", "ACAC", "AGAG", "ATTA", "TGAC", "AATT", "CGTT",
		"GTTC", "GGTA", "AGCA", "CATC"})
	actual = FrequentKmersWithMismatches(sequence, 4, 2)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}

	sequence, _ = CreateDNASequence("AATTAATTGGTAGGTAGGTA")
	expected, _ = DNASequences([]string{"GGTA"})
	actual = FrequentKmersWithMismatches(sequence, 4, 0)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected: ", expected, " Actual: ", actual)
	}
}

func TestRemoveDuplicates(t *testing.T) {
	sequences, err := DNASequences([]string{"ATCG", "ATCG", "ATGC", "ATGC", "GATC"})
	if err != nil {
		t.Error("didn't expect an error when creating a sequences")
	}

	expected, err := DNASequences([]string{"ATCG", "ATGC", "GATC"})
	if err != nil {
		t.Error("didn't expect an error when creating a sequences")
	}

	actual := RemoveDuplicates(sequences)
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	sequences, err = DNASequences([]string{"A", "G", "G"})
	if err != nil {
		t.Error("didn't expect an error when creating a sequences")
	}

	expected, err = DNASequences([]string{"A", "G"})
	if err != nil {
		t.Error("unexpected error when creating sequences")
	}

	actual = RemoveDuplicates(sequences)
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

func TestSequencesEqual(t *testing.T) {
	a, _ := DNASequences([]string{"A", "G"})
	b, _ := DNASequences([]string{"A", "G"})
	if !sequencesEqual(a, b) {
		t.Error("Expected a: ", a, " to equal b: ", b)
	}

	a, _ = DNASequences([]string{"G", "A"})
	b, _ = DNASequences([]string{"A", "G"})
	if !sequencesEqual(a, b) {
		t.Error("Expected a: ", a, " to equal b: ", b)
	}

	a, _ = DNASequences([]string{"A", "A"})
	b, _ = DNASequences([]string{"A", "A"})
	if !sequencesEqual(a, b) {
		t.Error("Expected a: ", a, " to equal b: ", b)
	}

	a, _ = DNASequences([]string{"A", "TA"})
	b, _ = DNASequences([]string{"A", "TA"})
	if !sequencesEqual(a, b) {
		t.Error("Expected a: ", a, " to equal b: ", b)
	}

	a, _ = DNASequences([]string{"A", "TA"})
	b, _ = DNASequences([]string{"A", "T"})
	if sequencesEqual(a, b) {
		t.Error("Expected a: ", a, " to not equal b: ", b)
	}
}

func TestSimilarPatterns(t *testing.T) {
	expected := []int{6, 7, 26, 27}
	genome, _ := CreateDNASequence("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT")
	pattern, _ := CreateDNASequence("ATTCTGGA")
	tolerance := 3
	actual := SimilarPatterns(genome, pattern, tolerance)

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	expected = []int{4, 5, 6, 7, 8, 11, 12, 13, 14, 15}
	genome, _ = CreateDNASequence("TTTTTTAAATTTTAAATTTTTT")
	pattern, _ = CreateDNASequence("AAA")
	tolerance = 2
	actual = SimilarPatterns(genome, pattern, tolerance)

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	expected = []int{0, 30, 66}
	genome, _ = CreateDNASequence("GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGCGCGCCGCGTCGATT")
	pattern, _ = CreateDNASequence("GAGCGCTGG")
	tolerance = 2
	actual = SimilarPatterns(genome, pattern, tolerance)

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	expected = []int{0, 1, 2, 3}
	genome, _ = CreateDNASequence("AAAAAA")
	pattern, _ = CreateDNASequence("TTT")
	tolerance = 3
	actual = SimilarPatterns(genome, pattern, tolerance)
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	expected = []int{0}
	genome, _ = CreateDNASequence("CCACCT")
	pattern, _ = CreateDNASequence("CCA")
	tolerance = 0
	actual = SimilarPatterns(genome, pattern, tolerance)
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	expected = []int{0, 1}
	genome, _ = CreateDNASequence("AAA")
	pattern, _ = CreateDNASequence("CA")
	tolerance = 1
	actual = SimilarPatterns(genome, pattern, tolerance)
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected ", expected, "Actual: ", actual)
	}
}

func TestPatternCount(t *testing.T) {
	genome, _ := CreateDNASequence("GCGCG")
	pattern, _ := CreateDNASequence("GCG")

	count := genome.PatternCount(pattern)
	if count != 2 {
		t.Error("Expected 2, got ", count)
	}
}

func TestPatternCountWithMismatches(t *testing.T) {
	genome, _ := CreateDNASequence("GCGCA")
	pattern, _ := CreateDNASequence("GCG")

	count := genome.PatternCountWithMismatches(pattern, 1)
	if count != 2 {
		t.Error("Expected 2, got ", count)
	}
}

// Checks to see if two sequences have the exact same dnaSequences but any order
func sequencesEqual(a []dnaSequence, b []dnaSequence) bool {
	if len(a) != len(b) {
		return false
	}
	for _, sequence := range a {
		var found bool = false
		for _, compareSequence := range b {
			if reflect.DeepEqual(sequence, compareSequence) {
				found = true
			}
		}
		if found == false {
			return false
		}
	}

	return true
}
