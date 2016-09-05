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
