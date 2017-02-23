from functools import reduce

def score(motifs):
    return map(reduce(lambda score, cur: score + cur, motifColumn), rotate(motifs))

def score(motifs):
    return map(reduce(
        lambda columnScore, nucleotide: columnScore.addToScore(nucleotide),
        rotate(motifs), MotifColumnScore()), rotate(motifs))


# Rotates the motifs 90 degrees
def rotate(motifs):
    tmp = [[]]
    for i in range(len(motifs[0])):
        for k in range(len(motifs)):
            tmp[k][i] = motifs[k][i]

class MotifColumnScore:
    score = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    def addToScore(nucleotide):
        if hasattr(score, nucleotide):
            score[nucleotide] += 1
        else:
            raise Exception("invalid nucleotide", nucleotide)

# find mininmum distance kmer in string
def minDistance(pattern, string):
    distances = []
    for i in len(string):
        currentMer = string[i:len(pattern)]
        distances[i] = distance(currentMer)

    return min(distances)


def medianString(k, patterns):
