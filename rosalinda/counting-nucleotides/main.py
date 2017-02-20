# dnaString is a string of nucleotides that will have each individual
# nucldeotide counted
def count(dnaString):
    nucleotides = {}
    for i in dnaString:
        if i in nucleotides:
            nucleotides[i] += 1
        else:
            nucleotides[i] = 1

    return nucleotides

dataset = open("dataset", "r").readline()
print(count(dataset))
