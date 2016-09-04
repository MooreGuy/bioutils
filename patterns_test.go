package bioutils

import "testing"

func TestPatternCount(t *testing.T) {
	var text []byte = []byte{'G', 'C', 'G', 'C', 'G'}
	var pattern string = "GCG"

	count := PatternCount(text, pattern)
	if count != 2 {
		t.Error("Expected 2, got ", count)
	}
}
