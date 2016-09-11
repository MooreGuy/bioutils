package bioutils

import (
	"bytes"
)

// returns the number of occurrences of pattern in a text
func PatternCount(text []byte, pattern string) int {
	var count int = 0
	var index int = 0
	for index <= len(text)-len(pattern) {
		if bytes.Equal(text[index:index+len(pattern)], []byte(pattern)) {
			count++
		}
		index++
	}

	return count
}

// returns the most frequent kmers in a text. Possible to have multiple
// most frequent
func MostFrequentKmers(text []byte, k int) []string {
	mostFrequentKmers := []string{}
	kmerPatternCounts := make([]int, len(text))
	maxCount := 0
	for i := 0; i <= len(text)-k; i++ {
		kmerPatternCounts[i] = PatternCount(text, string(text[i:i+k]))
		if kmerPatternCounts[i] > maxCount {
			maxCount = kmerPatternCounts[i]
		}
	}
	for i := 0; i <= len(text)-k; i++ {
		if kmerPatternCounts[i] == maxCount {
			mostFrequentKmers = append(mostFrequentKmers, string(text[i:i+k]))
		}
	}

	mostFrequentKmers = RemoveDuplicates(mostFrequentKmers)

	return mostFrequentKmers
}

func RemoveDuplicates(items []string) []string {
	for i := 0; i <= len(items)-1; i++ {
		for comparingIndex := i + 1; comparingIndex <= len(items)-1; comparingIndex++ {
			if items[i] == items[comparingIndex] {
				// delete i
				items = append(items[:comparingIndex], items[comparingIndex+1:]...)
				comparingIndex--

			}
		}
	}

	return items
}

// G to C skew
func Skew(genome string) []int {
	var skews []int = make([]int, len(genome)+1)
	var currentSkew int = 0

	skews[0] = currentSkew
	for i, nucleotide := range genome {
		if nucleotide == 'G' {
			currentSkew++
		} else if nucleotide == 'C' {
			currentSkew--
		}
		skews[i+1] = currentSkew
	}

	return skews
}

// Returns the lowest number and the indeces in the array of skews.
func FindLowest(skews []int) (lowest int, indexes []int) {
	lowest = skews[0]
	indexes = []int{0}
	for i := 1; i < len(skews); i++ {
		if lowest > skews[i] {
			lowest = skews[i]
			indexes = []int{i}
		} else if lowest == skews[i] {
			indexes = append(indexes, i)
		}
	}

	return lowest, indexes
}
