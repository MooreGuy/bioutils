import dnaSequence
import sys

def readSingle():
    string = open("dataset", "r").readline()
    return string.replace("\n", "")

def readMulti():
    dataset = []
    for i in iter(open("dataset", "r")):
        dataset.append(i.replace("\n", ""))

    return dataset

# Reads a single line from dataset and returns the resulting count
def count_nucleotides():
    count = dnaSequence.count(readSingle())
    print(count)

def transcribe_dna_to_rna():
   transcription = dnaSequence.transcribe(readSingle())
   print(transcription)

def reverse_complement():
    dnaString = readSingle()
    print(dnaSequence.complement(dnaString))

def distance():
    dnaStrings = readMulti()
    print(dnaSequence.distance(dnaStrings[0], dnaStrings[1]))

problems = {
    "count-nucleotides":     count_nucleotides,
    "transcribe-dna-to-rna": transcribe_dna_to_rna,
    "reverse-complement":    reverse_complement,
    "distance":              distance,
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
