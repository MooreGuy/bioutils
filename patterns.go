package bioutils

// Returns the lowest number and the indeces of the lowest in the array of
// numbers
func FindLowest(nums []int) (lowest int, indexes []int) {
	lowest = nums[0]
	indexes = []int{0}
	for i := 1; i < len(nums); i++ {
		if lowest > nums[i] {
			lowest = nums[i]
			indexes = []int{i}
		} else if lowest == nums[i] {
			indexes = append(indexes, i)
		}
	}

	return lowest, indexes
}

/*
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
*/

// I really don't like this, but since SimilarPatterns returns start indexes
// rather than the kmer itself, I have to do this.
/*
func RetrieveKmersFromIndexSlice(text []byte, startIndexes []int, k int) (kmers []string) {
	kmers = []string{}
	for _, startIndex := range startIndexes {
		kmers = append(kmers, string(text[startIndex:startIndex+k]))
	}

	return
}
*/

/*
func FindAllKmers(text []byte, k int) []string {
	for i := 0; i <= len(text)-k; i++ {
		kmers = append(kmers, string(text[i:i+k]))
	}
	return RemoveDuplicates(kmers)
}
*/
