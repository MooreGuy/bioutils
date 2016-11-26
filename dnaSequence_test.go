package bioutils

import (
	"log"
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

func testDiff(t *testing.T) {
	sequencesA, _ := DNASequences([]string{"AGT", "AGA"})
	sequencesB, _ := DNASequences([]string{"AGT", "AGC"})
	expectedDiffA, _ := DNASequences([]string{"AGA"})
	expectedDiffB, _ := DNASequences([]string{"AGC"})
	actualDiffA, actualDiffB := sequenceGroup(sequencesA).Diff(sequencesB)
	if !sequencesEqual(expectedDiffA, actualDiffA) {
		t.Error("Diff A doesn't match expected.")
	}
	if !sequencesEqual(expectedDiffB, actualDiffB) {
		t.Error("Diff B doesn't match expected.")
	}

	sequencesA, _ = DNASequences([]string{"AGT", "AGT"})
	sequencesB, _ = DNASequences([]string{"AGT", "AGT", "AGC"})
	expectedDiffA, _ = DNASequences([]string{})
	expectedDiffB, _ = DNASequences([]string{"AGC"})
	actualDiffA, actualDiffB = sequenceGroup(sequencesA).Diff(sequencesB)
	if !sequencesEqual(expectedDiffA, actualDiffA) {
		t.Error("Diff A doesn't match expected.")
	}
	if !sequencesEqual(expectedDiffB, actualDiffB) {
		t.Error("Diff B doesn't match expected.")
	}
}

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

	sequences, err = DNASequences([]string{"A", "G", "G", "G", "G"})
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

func TestRemove(t *testing.T) {
	group, _ := DNASequences([]string{"A", "G", "T"})
	removeIndex := 2
	expected, _ := DNASequences([]string{"A", "G"})
	actual := sequenceGroup(group).Remove(removeIndex)
	if !sequencesEqual(expected, actual) {
		t.Error("Remove didn't remove correctly.")
	}

	group, _ = DNASequences([]string{"A", "G", "T", "C"})
	removeIndex = 0
	expected, _ = DNASequences([]string{"G", "T", "C"})
	actual = sequenceGroup(group).Remove(removeIndex)
	if !sequencesEqual(expected, actual) {
		t.Error("Remove didn't remove correctly.")
	}

	group, _ = DNASequences([]string{"A", "G", "T", "C", "C", "C"})
	removeIndex = 3
	expected, _ = DNASequences([]string{"A", "G", "T", "C", "C"})
	actual = sequenceGroup(group).Remove(removeIndex)
	if !sequencesEqual(expected, actual) {
		t.Error("Remove didn't remove correctly.")
	}
}

func TestFindMotifs(t *testing.T) {
	sequence1, _ := CreateDNASequence("ATTTGGC")
	sequence2, _ := CreateDNASequence("TGCCTTA")
	sequence3, _ := CreateDNASequence("CGGTATC")
	sequence4, _ := CreateDNASequence("GAAAATT")
	sequences := []dnaSequence{sequence1, sequence2, sequence3, sequence4}
	expected, _ := DNASequences([]string{"ATA", "ATT", "GTT", "TTT"})
	actual := FindMotifs(sequences, 3, 1)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected motifs didn't match actual.")
	}

	sequences, _ = DNASequences([]string{"ACGT", "ACGT", "ACGT"})
	expected, _ = DNASequences([]string{"ACG", "CGT"})
	actual = FindMotifs(sequences, 3, 0)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected motifs didn't match actual.")
	}

	sequences, _ = DNASequences([]string{"AAAAA", "AAAAA", "AAAAA"})
	expected, _ = DNASequences([]string{"AAA", "AAC", "AAG", "AAT", "ACA", "AGA", "ATA", "CAA", "GAA", "TAA"})
	actual = FindMotifs(sequences, 3, 1)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected motifs didn't match actual.")
	}

	sequences, _ = DNASequences([]string{"AAAAA", "AAAAA", "AACAA"})
	expected, _ = DNASequences([]string{})
	actual = FindMotifs(sequences, 3, 0)
	if !sequencesEqual(expected, actual) {
		t.Error("Expected motifs didn't match actual.")
	}

	sequences, _ = DNASequences([]string{"AAAAA", "AAAAA", "AAAAA"})
	expected, _ = DNASequences([]string{"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"})
	actual = FindMotifs(sequences, 3, 3)
	if !sequencesEqual(expected, actual) {
		log.Println("actual", sequenceGroup(actual).Sort())
		log.Println("expected", sequenceGroup(expected).Sort())
		diffA, diffB := sequenceGroup(actual).Diff(actual)
		log.Println("diffA", diffA)
		log.Println("diffB", diffB)
		t.Error("Expected motifs didn't match actual.")
	}

}

func TestAllContain(t *testing.T) {
	sequence1, _ := CreateDNASequence("GATG")
	sequence2, _ := CreateDNASequence("GATG")
	sequences := []dnaSequence{sequence1, sequence2}
	toFind, _ := CreateDNASequence("AT")
	expected := true
	actual := AllContain(sequences, toFind, 0)
	if expected != actual {
		t.Error("Expected to find in all sequences, but didn't")
	}

	toFind, _ = CreateDNASequence("GT")
	expected = true
	actual = AllContain(sequences, toFind, 1)
	if expected != actual {
		t.Error("Expected to find in all sequences, but didn't")
	}

	toFind, _ = CreateDNASequence("CAT")
	expected = true
	actual = AllContain(sequences, toFind, 1)
	if expected != actual {
		t.Error("Expected to find in all sequences, but didn't")
	}

	toFind, _ = CreateDNASequence("CAT")
	expected = false
	actual = AllContain(sequences, toFind, 0)
	if expected != actual {
		t.Error("Expected to find in all sequences, but didn't")
	}

	toFind, _ = CreateDNASequence("CAG")
	expected = false
	actual = AllContain(sequences, toFind, 1)
	if expected != actual {
		t.Error("Expected to find in all sequences, but didn't")
	}
}

func TestContains(t *testing.T) {
	sequence, _ := CreateDNASequence("ATTTGGC")
	toFind, _ := CreateDNASequence("ATTT")
	expected := true
	actual := sequence.Contains(toFind, 0)
	if expected != actual {
		t.Error("Expected to find sequence, but didn't")
	}

	toFind, _ = CreateDNASequence("AGTT")
	expected = true
	actual = sequence.Contains(toFind, 1)
	if expected != actual {
		t.Error("Expected to find sequence, but didn't")
	}

	toFind, _ = CreateDNASequence("AGG")
	expected = false
	actual = sequence.Contains(toFind, 0)
	if expected != actual {
		t.Error("Didn't expect to find sequence, but did.")
	}

	toFind, _ = CreateDNASequence("AGGT")
	expected = false
	actual = sequence.Contains(toFind, 1)
	if expected != actual {
		t.Error("Didn't expect to find sequence, but did.")
	}

}

func testFindSequence(t *testing.T) {
	group, _ := DNASequences([]string{"AGT", "ACT", "GAC"})
	toFind, _ := CreateDNASequence("AGT")
	expected := 1
	actual := sequenceGroup(group).FindSequence(toFind)
	if expected != actual {
		t.Error("Expected index didn't match actual")
	}

	toFind, _ = CreateDNASequence("GCC")
	expected = -1
	actual = sequenceGroup(group).FindSequence(toFind)
	if expected != actual {
		t.Error("Expected index didn't match actual")
	}

	toFind, _ = CreateDNASequence("GAC")
	expected = 2
	actual = sequenceGroup(group).FindSequence(toFind)
	if expected != actual {
		t.Error("Expected index didn't match actual")
	}

}

func testEquals(t *testing.T) {
	sequenceA, _ := CreateDNASequence("AG")
	sequenceB, _ := CreateDNASequence("AG")
	expected := true
	actual := sequenceA.Equals(sequenceB)
	if expected != actual {
		t.Error("sequence A doesn't equal sequence B")
	}

	sequenceB, _ = CreateDNASequence("AC")
	expected = false
	actual = sequenceA.Equals(sequenceB)
	if expected != actual {
		t.Error("sequence A doesn't equal sequence B")
	}

	sequenceB, _ = CreateDNASequence("AGC")
	expected = false
	actual = sequenceA.Equals(sequenceB)
	if expected != actual {
		t.Error("sequence A doesn't equal sequence B")
	}

}

// Checks to see if two sequence slices have the exact same dnaSequences but
// any order
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
