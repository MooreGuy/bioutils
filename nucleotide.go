package bioutils

// Represents a single monomer, and is identified by the base letter.
type nucleotide rune
type padNucleotide nucleotide

const padLetter = 'F'

func PadNucleotide() padNucleotide {
	return padLetter
}

func (n nucleotide) ValidNucleotide() bool {
	_, ok := GetValidNucleotidesMap()[rune(n)]
	return ok
}

// Returns all valid nucleotides within a map for quick searching.
func GetValidNucleotidesMap() map[rune]bool {
	return map[rune]bool{
		'A': true,
		'T': true,
		'C': true,
		'G': true,
	}
}

// Returns all nucleotides in an array
func GetValidNucleotides() [4]nucleotide {
	return [4]nucleotide{'A', 'T', 'C', 'G'}
}

// Returns all valid nucleotides that aren't this one.
func (n nucleotide) OtherNucleotides() (otherNucleotides []nucleotide) {
	index := 0
	for currentNucleotide, _ := range GetValidNucleotidesMap() {
		if currentNucleotide != rune(n) {
			otherNucleotides[index] = nucleotide(currentNucleotide)
			index++
		}
	}
	return
}
