def transcribe(dnaString):
    transcription = ""
    for i in dnaString:
        if i == "T":
            transcription += "U"
        else:
            transcription += i

    print(transcription)
    return transcription

dataset = open("dataset", "r")
output = open("output", "w")
output.write(transcribe(dataset.readline()))
