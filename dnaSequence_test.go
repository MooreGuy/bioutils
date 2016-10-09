package bioutils

import (
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
