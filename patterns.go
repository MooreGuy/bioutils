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

// Finds patterns similar within a tolerance and returns their starting indexes
func SimilarPatterns(text []byte, pattern string, tolerance int) []int {
	patternIndexes := []int{}
	for i := 0; i < len(text); i++ {
		endIndex := min(i+len(pattern), len(text)-1)
		compareString := string(text[i:endIndex])
		if len(compareString) < len(pattern) && i < len(pattern) {
			compareString = LeftPad(compareString, len(pattern)-len(compareString))
		}
		if HammingDistance(compareString, pattern) <= tolerance {
			patternIndexes = append(patternIndexes, i)
		}
	}

	return patternIndexes
}

func LeftPad(toPad string, padLength int) string {
	padding := make([]byte, 0, padLength)
	for i := 0; i < padLength; i++ {
		padding = append(padding, 'F')
	}

	return string(append(padding, []byte(toPad)...))
}

func min(x, y int) int {
	if x <= y {
		return x
	}
	return y
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

func MostFrequentKmersWithMismatch(text []byte, k int, tolerance int) []string {
	allKmers := FindAllKmers(text, k)
	maxCountKmers := []string{}
	maxCount := 0
	for _, kmer := range allKmers {
		similarPatternIndexes := SimilarPatterns(text, kmer, tolerance)
		if len(similarPatternIndexes) > maxCount {
			maxCountKmers = []string{kmer}
			// Retrieve all similar kmers.
			maxCountKmers = append(maxCountKmers,
				RetrieveKmersFromIndexSlice(text, similarPatternIndexes, k)...)
			maxCount = len(similarPatternIndexes)
		} else if len(similarPatternIndexes) == maxCount {
			maxCountKmers = append(maxCountKmers, kmer)
		}
	}

	return RemoveDuplicates(maxCountKmers)
}

// I really don't like this, but since SimilarPatterns returns start indexes
// rather than the kmer itself, I have to do this.
func RetrieveKmersFromIndexSlice(text []byte, startIndexes []int, k int) (kmers []string) {
	kmers = []string{}
	for _, startIndex := range startIndexes {
		kmers = append(kmers, string(text[startIndex:startIndex+k]))
	}

	return
}

func FindAllKmers(text []byte, k int) (kmers []string) {
	for i := 0; i <= len(text)-k; i++ {
		kmers = append(kmers, string(text[i:i+k]))
	}
	kmers = RemoveDuplicates(kmers)
	return
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

// The `5 to `3 skew or G to C skew
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

func HammingDistance(p string, q string) int {
	mismatches := 0
	var shorterString string
	var longerString string
	if len(p) > len(q) {
		shorterString = q
		longerString = p
	} else {
		shorterString = p
		longerString = q
	}

	for i, value := range shorterString {
		if value != rune(longerString[i]) {
			mismatches++
		}
	}

	mismatches += len(longerString) - len(shorterString)
	return mismatches
}
