import dnaSequence
import sys

def readSingleNucleotide():
    string = open("dataset", "r").readline()
    return string.replace("\n", "")

# Reads a single line from dataset and returns the resulting count
def count_nucleotides():
    count = dnaSequence.count(readSingleNucleotide())
    print(count)

def transcribe_dna_to_rna():
   transcription = dnaSequence.transcribe(readSingleNucleotide())
   print(transcription)

def reverse_complement():
    dnaString = readSingleNucleotide()
    print(dnaSequence.complement(dnaString))

problems = {
    "count-nucleotides": count_nucleotides,
    "transcribe-dna-to-rna": transcribe_dna_to_rna,
    "reverse-complement": reverse_complement
}
print(count_nucleotides())

print("Select which problem to run:")
print("============================")
for key in problems:
    print(key)

prob = raw_input("Enter the desired problem: ")
if prob in problems:
    problems[prob]()
else:
    print("Couldn't find function for problem: " + prob)