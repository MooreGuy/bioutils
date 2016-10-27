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
