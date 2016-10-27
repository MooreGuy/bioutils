package bioutils

import (
	"errors"
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

	mostFrequentKmers.sequences = RemoveDuplicates(mostFrequentKmers.sequences)

	return mostFrequentKmers
}
