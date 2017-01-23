package bioutils

import (
	"reflect"
	"testing"
)

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
