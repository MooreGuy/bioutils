package bioutils

import (
	"reflect"
	"testing"
)

func TestPatternCount(t *testing.T) {
	var text []byte = []byte{'G', 'C', 'G', 'C', 'G'}
	var pattern string = "GCG"

	count := PatternCount(text, pattern)
	if count != 2 {
		t.Error("Expected 2, got ", count)
	}
}

func TestMostFrequentKmers(t *testing.T) {
	var text []byte = []byte("ACGTTGCATGTCGCATGATGCATGAGAGCT")
	var k int = 4

	frequentWords := MostFrequentKmers(text, k)
	expected := []string{"GCAT", "CATG"}
	if !reflect.DeepEqual(frequentWords, expected) {
		t.Error("Expected: ", expected, "Actual: ", frequentWords)
	}
}

func TestRemoveDuplicates(t *testing.T) {
	items := []string{"abba", "Abba", "abba", "Abba", "gabba"}
	expected := []string{"abba", "Abba", "gabba"}
	actual := RemoveDuplicates(items)
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

func TestSkew(t *testing.T) {
	expected := []int{0, -1, -1, -1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 1, 0, 0, 0, 0, -1, 0, -1, -2}
	skew := Skew("CATGGGCATCGGCCATACGCC")
	if !reflect.DeepEqual(expected, skew) {
		t.Error("Expected: ", expected, "Actual: ", skew)
	}

	skew = Skew("GAGCCACCGCGATA")
	expected = []int{0, 1, 1, 2, 1, 0, 0, -1, -2, -1, -2, -1, -1, -1, -1}

	if !reflect.DeepEqual(expected, skew) {
		t.Error("Expected: ", expected, "Actual: ", skew)
	}
}

func TestFindLowest(t *testing.T) {
	expected := []int{11, 24}
	_, actual := FindLowest(Skew("TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"))

	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

func TestHammingDistance(t *testing.T) {
	expected := 3
	actual := HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC")
	if expected != actual {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

func TestLeftPad(t *testing.T) {
	expected := "FFFABC"
	actual := LeftPad("ABC", 3)
	if expected != actual {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

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

func TestFindAllKmers(t *testing.T) {
	text := []byte("GACTGACT")
	k := 3
	actual := FindAllKmers(text, k)
	expected := []string{"GAC", "ACT", "CTG", "TGA"}
	if !reflect.DeepEqual(expected, actual) {
		t.Error("Expected: ", expected, "Actual: ", actual)
	}
}

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
	expected = MostFrequentKmersWithMismatch(text, k, tolerance)
	actual = []string{"ATGT", "GATG", "ATGC"}
	if !reflect.DeepEqual(actual, expected) {
		t.Error("Exepcted: ", expected, "Actual: ", actual)
	}
}
