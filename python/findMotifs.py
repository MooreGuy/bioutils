import dnaSequence

dataset = open("dataset", "r")
k = 0
d = 0
sequences = []
for i, line in enumerate(dataset):
    currentLine = line.strip()
    if i == 0:
        kandd = currentLine.split(" ")
        k = int(kandd[0])
        d = int(kandd[1])
    else:
        sequences.append(currentLine)

dnaSequence.output(dnaSequence.findMotifs(sequences,k,d))
