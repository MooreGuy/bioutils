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
