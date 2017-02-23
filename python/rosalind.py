import dnaSequence
import sys

def readSingleNucleotide():
    return open("dataset", "r").readline()

# Reads a single line from dataset and returns the resulting count
def count_nucleotides():
    count = dnaSequence.count(readSingleNucleotide())
    print(count)

def transcribe_dna_to_rna():
   transcription = dnaSequence.transcribe(readSingleNucleotide())
   print(transcription)

problems = {
    "count-nucleotides": count_nucleotides,
    "transcribe-dna-to-rna": transcribe_dna_to_rna,
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
