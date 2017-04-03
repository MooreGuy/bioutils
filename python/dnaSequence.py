def kmers(dnaSequence, k):
    kmers = []
    for i in range(len(dnaSequence)-k+1):
        kmers.append(dnaSequence[i:i+k])
    return kmers

def neighbors(dnaSequence, d):
    if d == 0:
        return [dnaSequence]

    if (len(dnaSequence) == 1):
        return ["A", "C", "G", "T"]

    neighborhood = []

    suffix = dnaSequence[1:]
    suffixNeighbors = neighbors(suffix, d)
    for suffixNeighbor in suffixNeighbors:
        # If the hamming distance allows a mismatch, add all possible mismatches
        if (distance(suffix, suffixNeighbor) < d):
            for nucleotide in ["A", "C", "G", "T"]:
                newSequence = nucleotide + suffixNeighbor
                neighborhood.append(newSequence)
        else:
            firstNucleotide = dnaSequence[0]
            newSequence = firstNucleotide + suffixNeighbor
            neighborhood.append(newSequence)

    return neighborhood

def distance(sequenceA, sequenceB):
    mismatches = 0
    if len(sequenceA) > len(sequenceB):
        longerSequence = sequenceA
        shorterSequence = sequenceB
    else:
        longerSequence = sequenceB
        shorterSequence = sequenceA
    for i in range(len(shorterSequence)):
        if longerSequence[i] != shorterSequence[i]:
            mismatches += 1
    mismatches += len(longerSequence) - len(shorterSequence)
    return mismatches

def findMotifs(sequences, k, d):
    motifs = []
    mers = []
    for sequence in sequences:
        mers.extend(kmers(sequence, k))
    for mer in mers:
        for neighbor in neighbors(mer, d):
            if isMotif(sequences, neighbor, d):
                motifs.append(neighbor)

    return set(motifs)

def isMotif(sequences, pattern, d):
    for sequence in sequences:
        if contains(sequence, pattern, d) == False:
            return False

    return True

def contains(sequence, pattern, d):
    for i in range(0, (len(sequence) - len(pattern)) + 1):
        if distance(sequence[i:i+len(pattern)], pattern) <= d:
            return True
    return False

def output(sequences):
    for sequence in sequences:
        print(sequence + " ")

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

# Tates a DNA string and transcribes it to RNA
def transcribe(dnaString):
    transcription = ""
    for i in dnaString:
        if i == "T":
            transcription += "U"
        else:
            transcription += i

    print(transcription)
    return transcription

# Given a DNA string gives the reverse complement
# ex: ACG -> CGT
def complement(dnaString):
    complements = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }

    revComplement = map(lambda nuc: {complements[nuc]}, dnaString[::-1])
    # map() returns a list of sets? Maybe my doing, but this will fix it
    return ''.join(i.pop() for i in revComplement)

# Given a string and a potential substring, returns all
def subs(string, subs):
    subsLocations = []
    for i in range(len(string)):
        if subs == string[i:i + len(subs)]:
            subsLocations.append(i+1)

    return subsLocations
