package bioutils

import (
	"errors"
	"reflect"
)

// Represents a sequence of nucleotides
type dnaSequence []nucleotide

// Takes in a string representation of the nucleotide bases and returns a
// completely checked dna sequence
func CreateDNASequence(uncheckedSequence string) (dnaSequence, error) {
	newSequence := []nucleotide{}
	for _, currentBaseLetter := range uncheckedSequence {
		curNucleotide := nucleotide(currentBaseLetter)
		if !curNucleotide.ValidNucleotide() {
			return newSequence, errors.New("Invalid base " + string(currentBaseLetter))
		}
		newSequence = append(newSequence, curNucleotide)
	}

	return newSequence, nil
}

func DNASequences(uncheckedSequences []string) ([]dnaSequence, error) {
	sequences := []dnaSequence{}
	for _, unchecked := range uncheckedSequences {
		checked, err := CreateDNASequence(unchecked)
		if err != nil {
			return sequences, err
		}

		sequences = append(sequences, checked)
	}

	return sequences, nil
}

func DNASequencesFromNucleotides(nucleotides []nucleotide) []dnaSequence {
	sequence := []dnaSequence{}
	for _, nucleotide := range nucleotides {
		sequence = append(sequence, dnaSequence{nucleotide})
	}
	return sequence
}

// The hamming distance is the number of differences of two dnaSequences
func (sequence dnaSequence) HammingDistance(compareSequence dnaSequence) int {
	mismatches := 0
	var shorterSequence dnaSequence
	var longerSequence dnaSequence
	if len(sequence) > len(compareSequence) {
		shorterSequence = compareSequence
		longerSequence = sequence
	} else {
		shorterSequence = sequence
		longerSequence = compareSequence
	}

	for i, value := range shorterSequence {
		if value != longerSequence[i] {
			mismatches++
		}
	}

	mismatches += len(longerSequence) - len(shorterSequence)
	return mismatches
}

// Finds patterns similar within a tolerance and returns their starting indexes
func SimilarPatterns(genome dnaSequence, pattern dnaSequence, tolerance int) []int {
	patternIndexes := []int{}
	for i := 0; i < len(genome); i++ {
		endIndex := min(i+len(pattern), len(genome)-1)
		comparePattern := genome[i:endIndex]
		if len(comparePattern) < len(pattern) && i < len(pattern) {
			comparePattern = LeftPad(comparePattern, len(pattern)-len(comparePattern))
		}
		if comparePattern.HammingDistance(pattern) <= tolerance {
			patternIndexes = append(patternIndexes, i)
		}
	}

	return patternIndexes
}

// Left pad for Hamming distance, pads with bad nucleotides
func LeftPad(toPad dnaSequence, padLength int) dnaSequence {
	padding := make([]nucleotide, 0, padLength)
	for i := 0; i < padLength; i++ {
		padding = append(padding, nucleotide(PadNucleotide()))
	}

	return append(padding, toPad...)
}

func min(x, y int) int {
	if x <= y {
		return x
	}
	return y
}

func (genome dnaSequence) PatternCount(pattern dnaSequence) int {
	var count int = 0
	var index int = 0
	for index <= len(genome)-len(pattern) {
		if reflect.DeepEqual(genome[index:index+len(pattern)], pattern) {
			count++
		}
		index++
	}

	return count
}

func (sequence dnaSequence) String() string {
	byteSequence := []byte{}
	for _, nucleotide := range sequence {
		byteSequence = append(byteSequence, byte(nucleotide))
	}
	return string(byteSequence)
}

// The `5 to `3 skew or G to C skew
func Skew(genome dnaSequence) []int {
	var skews []int = make([]int, len(genome)+1)
	var currentSkew int = 0

	skews[0] = currentSkew
	for i, curNucleotide := range genome {
		if curNucleotide == 'G' {
			currentSkew++
		} else if curNucleotide == 'C' {
			currentSkew--
		}
		skews[i+1] = currentSkew
	}

	return skews
}

// A neighbor  are kmers that are within a Hamming Distance away from the given
// pattern.
func GenerateNeighbors(pattern dnaSequence, distance int) (neighborhood []dnaSequence) {
	neighborhood = []dnaSequence{}
	if distance == 0 {
		neighborhood = []dnaSequence{pattern}
		return
	}
	if len(pattern) == 1 {
		neighborhood = DNASequencesFromNucleotides(GetValidNucleotidesSlice())
		return
	}

	suffix := pattern[1:]
	suffixNeighbors := GenerateNeighbors(suffix, distance)
	for _, sequence := range suffixNeighbors {
		// If the hamming distance allows a mismatch, add all possible mismatches
		if suffix.HammingDistance(sequence) < distance {
			for _, nucl := range GetValidNucleotides() {
				newSequence := append(dnaSequence{nucl}, sequence...)
				neighborhood = append(neighborhood, newSequence)
			}
		} else {
			firstNucleotideSequence := dnaSequence{pattern[0]}
			newSequence := append(firstNucleotideSequence, sequence...)
			neighborhood = append(neighborhood, newSequence)
		}
	}

	return
}
