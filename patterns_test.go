package bioutils

import (
	"reflect"
	"testing"
)

/*
func TestPatternCount(t *testing.T) {
	var text []byte = []byte{'G', 'C', 'G', 'C', 'G'}
	var pattern string = "GCG"

	sequence1, err := CreateDNASequence("GCGCG")
	if err != nil {
		t.Error("didn't expect an error when creating a sequence")
	}

	sequence2, err := CreateDNASequence("GCG")
	if err != nil {
		t.Error("didn't expect an error when creating a sequence")
	}

	count := sequence1.PatternCount(sequence2)
	if count != 2 {
		t.Error("Expected 2, got ", count)
	}
}
*/

func TestMostFrequentKmers(t *testing.T) {
	genome, err := CreateDNASequence("ACGTTGCATGTCGCATGATGCATGAGAGCT")
	if err != nil {
		t.Error("didn't expect an error when creating a sequence")
	}
	var k int = 4

	frequentWords := MostFrequentKmers(genome, k)
	sequence1, err := CreateDNASequence("GCAT")
	if err != nil {
		t.Error("didn't expect an error when creating a sequence")
	}

	sequence2, err := CreateDNASequence("CATG")
	if err != nil {
		t.Error("didn't expect an error when creating a sequence")
	}

	expectedSequences := []dnaSequence{sequence1, sequence2}

	if !reflect.DeepEqual(expectedSequences, frequentWords.Sequences()) {
		t.Error("Expected: ", expectedSequences, "Actual: ", frequentWords.Sequences())
	}
}

func TestRemoveDuplicates(t *testing.T) {
	sequences, err := DNASequences([]string{"ATCG", "ATCG", "ATGC", "ATGC", "GATC"})
	if err != nil {
		t.Error("didn't expect an error when creating a sequences")
	}
	kmerWithDuplicates, err := CreateKmer(sequences)
	if err != nil {
		t.Error("didn't expect an error when creating kmer")
	}

	expected, err := DNASequences([]string{"ATCG", "ATGC", "GATC"})
	if err != nil {
		t.Error("didn't expect an error when creating a sequences")
	}

	actual := kmerWithDuplicates.RemoveDuplicates()
	if !reflect.DeepEqual(expected, actual.Sequences()) {
		t.Error("Expected: ", expected, "Actual: ", actual.Sequences())
	}
}

func TestSkew(t *testing.T) {
	expected := []int{0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2}
	sequence, err := CreateDNASequence("CATGGGCATCGGCCATACGCC")
	if err != nil {
		t.Error("didn't expect an error when creating a sequence")
	}
	skew := Skew(sequence)
	if !reflect.DeepEqual(expected, skew) {
		t.Error("Expected: ", expected, "Actual: ", skew)
	}

	sequence, err = CreateDNASequence("GAGCCACCGCGATA")
	if err != nil {
		t.Error("didn't expect an error when creating a sequence")
	}
	skew = Skew(sequence)
	expected = []int{0, 1, 1, 2, 1, 0, 0, -1, -2, -1, -2, -1, -1, -1, -1}

	if !reflect.DeepEqual(expected, skew) {
		t.Error("Expected: ", expected, "Actual: ", skew)
	}
}

func TestFindLowest(t *testing.T) {
	expected := []int{11, 24}
	sequence, err := CreateDNASequence("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT")
	if err != nil {
		t.Error(err.Error())
	}
	_, actual := FindLowest(Skew(sequence))

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

/*
func TestHammingDistance(t *testing.T) {
	expected := 3
	actual := HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC")
	if expected != actual {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}
*/

/*
func TestLeftPad(t *testing.T) {
	expected := "FFFABC"
	actual := LeftPad("ABC", 3)
	if expected != actual {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}
*/

/*
func TestSimilarPatterns(t *testing.T) {
	expected := []int{6, 7, 26, 27}
	text := []byte("CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT")
	pattern := "ATTCTGGA"
	tolerance := 3
	actual := SimilarPatterns(text, pattern, tolerance)

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	expected = []int{4, 5, 6, 7, 8, 11, 12, 13, 14, 15}
	text = []byte("TTTTTTAAATTTTAAATTTTTT")
	pattern = "AAA"
	tolerance = 2
	actual = SimilarPatterns(text, pattern, tolerance)

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}

	expected = []int{0, 30, 66}
	text = []byte("GAGCGCTGGGTTAACTCGCTACTTCCCGACGAGCGCTGTGGCGCAAATTGGCGATGAAACTGCAGAGAGAACTGGTCATCCAACTGAATTCTCCCCGCTATCGCATTTTGATGCGCGCCGCGTCGATT")
	pattern = "GAGCGCTGG"
	tolerance = 2
	actual = SimilarPatterns(text, pattern, tolerance)

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}
*/

/*
func TestFindAllKmers(t *testing.T) {
	text := []byte("GACTGACT")
	k := 3
	actual := FindAllKmers(text, k)
	expected := []string{"GAC", "ACT", "CTG", "TGA"}
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}
*/

/*
func TestMostFrequentKmersWithMismatch(t *testing.T) {
	text := []byte("GACTGACT")
	k := 3
	tolerance := 1
	actual := MostFrequentKmersWithMismatch(text, k, tolerance)
	expected := []string{"GAC", "ACT"}
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Exepcted: ", expected, "Actual: ", actual)
	}

	text = []byte("ACGTTGCATGTCGCATGATGCATGAGAGCT")
	k = 4
	tolerance = 1
	actual = MostFrequentKmersWithMismatch(text, k, tolerance)
	expected = []string{"ATGT", "GATG", "ATGC"}
	if !reflect.DeepEqual(actual, expected) {
		t.Error("Exepcted: ", expected, "Actual: ", actual)
	}

	text = []byte("AAAAAAAAAA")
	k = 2
	tolerance = 1
	actual = MostFrequentKmersWithMismatch(text, k, tolerance)
	expected = []string{"AA", "AC", "AG", "CA", "AT", "GA", "TA"}
	if !reflect.DeepEqual(actual, expected) {
		t.Error("Exepcted: ", expected, "Actual: ", actual)
	}
}
*/
