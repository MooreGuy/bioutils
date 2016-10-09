package bioutils

import (
	"errors"
	"reflect"
)

type kmer struct {
	k         int
	sequences []dnaSequence
}

func (thisKmer kmer) Sequences() []dnaSequence {
	return thisKmer.sequences
}

func CreateKmer(sequences []dnaSequence) (kmer, error) {
	if len(sequences) < 1 {
		return kmer{}, errors.New("Cannot create kmers from a null set")
	}

	expectedRange := len(sequences[0])
	for _, currentSequence := range sequences {
		if len(currentSequence) != expectedRange {
			return kmer{}, errors.New("All sequences need to be of length k")
		}
	}

	return kmer{sequences: sequences, k: expectedRange}, nil
}

// A neighbor  are kmers that are within a Hamming Distance away from the given
// pattern.
func GenerateNeighbors(pattern dnaSequence, distance int) (neighborhoodKmer kmer) {
	sequences := []dnaSequence{}
	sequences = append(sequences, pattern)
	neighborhoodKmer = kmer{k: len(pattern), sequences: sequences}
	if distance == 0 {
		return
	}
	if len(pattern) == 1 {
		sequences = append(sequences, pattern[0].OtherNucleotides())
		return
	}

	suffixNeighbors := GenerateNeighbors(pattern[1:], distance)
	for _, sequence := range suffixNeighbors.sequences {
		suffix := pattern[1:]
		if suffix.HammingDistance(sequence) < distance {
			for _, nucl := range GetValidNucleotides() {
				newSequence := append(dnaSequence{nucl}, sequence...)
				sequences = append(sequences, newSequence)
			}
		}
	}
	return
}

// returns the most frequent kmers in a genome. Possible to have multiple
// most frequent
func MostFrequentKmers(genome dnaSequence, k int) (mostFrequentKmers kmer) {
	mostFrequentKmers = kmer{sequences: []dnaSequence{}, k: k}
	kmerPatternCounts := make([]int, len(genome))
	maxCount := 0
	for i := 0; i <= len(genome)-k; i++ {
		kmerPatternCounts[i] = genome.PatternCount(genome[i : i+k])
		if kmerPatternCounts[i] > maxCount {
			maxCount = kmerPatternCounts[i]
		}
	}
	for i := 0; i <= len(genome)-k; i++ {
		if kmerPatternCounts[i] == maxCount {
			mostFrequentKmers.sequences = append(mostFrequentKmers.sequences, genome[i:i+k])
		}
	}

	mostFrequentKmers = mostFrequentKmers.RemoveDuplicates()

	return mostFrequentKmers
}

func (merDuplicates kmer) RemoveDuplicates() kmer {
	for i := 0; i <= len(merDuplicates.sequences)-1; i++ {
		for comparingIndex := i + 1; comparingIndex <= len(merDuplicates.sequences)-1; comparingIndex++ {
			currentMer := merDuplicates.sequences[i]
			compareMer := merDuplicates.sequences[comparingIndex]
			if reflect.DeepEqual(currentMer, compareMer) {
				// delete i
				merDuplicates.sequences = append(merDuplicates.sequences[:comparingIndex],
					merDuplicates.sequences[comparingIndex+1:]...)
				comparingIndex--
			}
		}
	}

	return merDuplicates
}
