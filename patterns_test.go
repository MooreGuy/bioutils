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
