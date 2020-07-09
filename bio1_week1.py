# This script contains all algorithms needed for bioinformatics1


import numpy as np


def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count + 1
    return count


def FrequentWords(Text, k):

    """
    Finds the most common pattern in a length of text
    """

    FrequentPatterns = []
    Count = dict()

    for i in range(len(Text)-k+1):
        # This gets the pattern
        Pattern = Text[i:i+k]
        # This stores the amount of times that pattern exists in dict
        Count[Pattern] = PatternCount(Text, Pattern)

    # Store max val
    maxCount = max(Count.values())

    # Loop through
    for item in Count.items():
        if item[1] == maxCount:
            FrequentPatterns.append(item[0])

    return FrequentPatterns


def ReverseCompliment(Text):

    revcomp = ""

    for c in Text:
        if c == "A":
            c = "T"
        elif c == "C":
            c = "G"
        elif c == "G":
            c = "C"
        elif c == "T":
            c = "A"

        revcomp = c + revcomp
    return revcomp


def PatternMatch(Text, Pattern):
    indicies = []
    str_format = ""
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            indicies.append(i)

    for i in indicies:
        str_format = str_format + str(i) + " "
    print(str_format)


def ClumpFindKak(Text, k, L, t):
    '''
    In a window of length L within the Text we will find a k-mer
    that occures t times THIS WORKS BUT ITS KAK SLOW

    '''
    # Test dataset output AAACCAGGTGG

    FrequentPatterns = []
    Count = dict()

    for i in range(len(Text)-L+1):

        for j in range(len(Text[i:i+L])-k+1):
            # This gets the pattern
            Pattern = Text[j:j+k]
            # This stores the amount of times that pattern exists in dict
            Count[Pattern] = PatternCount(Text, Pattern)

        # Store max val
        maxCount = max(Count.values())

        # Loop through
        for item in Count.items():
            if item[1] >= t:
                if item[0] not in FrequentPatterns:
                    FrequentPatterns.append(item[0])

    return FrequentPatterns


def ClumpFindFast(Text, k, L, t):
    """
    call pattern count once then look at the first and last genome
    use dictionaries
    """
    print("DOne")


def SymbolToNumber(symbol):
    if symbol == "A":
        return 0
    if symbol == "C":
        return 1
    if symbol == "G":
        return 2
    if symbol == "T":
        return 3


def NumberToSymbol(index):
    if index == 0:
        return "A"
    if index == 1:
        return "C"
    if index == 2:
        return "G"
    if index == 3:
        return "T"


def PatternToNumber(Pattern):
    if len(Pattern) == 0:
        return 0

    symbol = Pattern[len(Pattern)-1]
    prefix = Pattern[0:len(Pattern)-1]

    return 4 * PatternToNumber(prefix) + SymbolToNumber(symbol)


def NumberToPattern(index, k):

    if k == 1:
        return NumberToSymbol(index)

    prefixIndex = index // 4  # This gives the quotient
    r = index % 4  # This gives remainder
    symbol = NumberToSymbol(r)

    PrefixPattern = NumberToPattern(prefixIndex, k - 1)

    return PrefixPattern + symbol


def main():
    # read the data
    # dat = np.loadtxt("test_dataset.txt", dtype=("str"))
    # First challenge
    # print(PatternCount(dat[0], dat[1]))

    # Second challenge
    # print(FrequentWords(dat[0], int(dat[1])))

    # Third challenge
    # print(ReverseCompliment(str(dat)))

    # Forth Challange
    # PatternMatch(dat[1], dat[0])

    # Honors challenge
    # PatternMatch(str(dat), "CTTGATCAT")

    # Fifth challenge
    # with open('test_dataset.txt', 'r') as myfile:
    #     dat = myfile.read().replace('\n', ' ')
    # myfile.close()

    # dat = dat.split()

    # print(ClumpFind(dat[0], int(dat[1]), int(dat[2]), int(dat[3])))
    print(ClumpFind("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA", 5, 50, 4))


if __name__ == "__main__":
    main()