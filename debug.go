package bioutils

import (
	"fmt"
)

func DisplayDiff(groupA sequenceGroup, groupB sequenceGroup) {
	diffA, diffB := groupA.Diff(groupB)
	fmt.Print("Different in A: ")
	fmt.Println(diffA)
	fmt.Print("Different in B: ")
	fmt.Println(diffB)
}
