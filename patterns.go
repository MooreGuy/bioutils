package bioutils

import (
	"bytes"
)

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
