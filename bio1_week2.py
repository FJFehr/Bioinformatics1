# This script contains all algorithms needed for bioinformatics1


import numpy as np

from bio1_week1 import PatternCount, ReverseCompliment


def Skew(Genome):

    skew_dict = {1: 0}
    skew_val = 0

    for i in range(len(Genome)):
        # create an updating skew value
        if Genome[i] == "G":
            skew_val = skew_val + 1
        elif Genome[i] == "C":
            skew_val = skew_val - 1

        # Use position as key and skew value as value
        skew_dict[i+1] = skew_val

    # find min value
    min_val = min(skew_dict.values())

    # find all positions of that min value
    min_val_positions = []
    for item in skew_dict.items():
        if item[1] == min_val:
            min_val_positions.append(item[0])

    return min_val_positions


def HammingDistance(p, q):

    ham_dist = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            ham_dist += 1
    return ham_dist


def ApproxPatternMatch(Text, Pattern, k):
    indicies = []
    for i in range(len(Text)-len(Pattern)+1):

        ham_dist = HammingDistance(Text[i:i+len(Pattern)], Pattern)

        if ham_dist <= k:
            indicies.append(i)

    return indicies


def PatternCountD(Text, Pattern, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):

        ham_dist = HammingDistance(Text[i:i+len(Pattern)], Pattern)

        if ham_dist <= d :
            count += 1
    return count


def ImmediateNeighbours(Pattern):

    Neighbourhood = []

    for i in range(len(Pattern)):
        symbol = Pattern[i]
        allsymbols = ["A", "C", "T", "G"]

        for j in allsymbols:
            if j != symbol:
                # insert the different j into the string

                Neighbour = Pattern[0:i] + j + Pattern[i+1:len(Pattern)]


                Neighbourhood.append(Neighbour)

    Neighbourhood.append(Pattern)

    return Neighbourhood


def Neighbours(Pattern, d):

    if d == 0:
        return Pattern

    elif len(Pattern) == 1:
        return ["A", "C", "G", "T"]

    else:

        neighbourhood = []

        nucleotide = ["A", "C", "G", "T"]

        # print("suff pattern looks like ", Pattern[1:len(Pattern)])

        suffixNeighbours = Neighbours(Pattern[1:len(Pattern)], d)



        for i in suffixNeighbours:



            if HammingDistance(Pattern[1:len(Pattern)], i) < d:

                for j in nucleotide:
                    neighbourhood.append(j + i)

            else:
                neighbourhood.append(Pattern[0] + i)

        return neighbourhood


def FrequentWordsWithMismatchRC(Text, k, d):

    """

    """
    FrequentPatterns = []
    Count = dict()

    for i in range(len(Text)-k+1):

        Pattern = Text[i:i+k]
        temp_neighbourhood = Neighbours(Pattern, d)

        for neighbour in temp_neighbourhood:
            if neighbour in Count:
                Count[neighbour] += 1
            else:
                Count[neighbour] = 1

        #Adding reverse compliment
        for neighbour in temp_neighbourhood:
            rc_neighbour = ReverseCompliment(neighbour)
            if rc_neighbour in Count:
                Count[rc_neighbour] += 1
            else:
                Count[rc_neighbour] = 1




    # Store max val
    maxCount = max(Count.values())

    # Loop through
    for item in Count.items():
        if item[1] == maxCount:
            FrequentPatterns.append(item[0])

    return FrequentPatterns


def FrequentWordsWithMismatch(Text, k, d):

    """

    """
    FrequentPatterns = []
    Count = dict()

    for i in range(len(Text)-k+1):

        Pattern = Text[i:i+k]
        temp_neighbourhood = Neighbours(Pattern, d)

        for neighbour in temp_neighbourhood:
            if neighbour in Count:
                Count[neighbour] += 1
            else:
                Count[neighbour] = 1

    # Store max val
    maxCount = max(Count.values())

    # Loop through
    for item in Count.items():
        if item[1] == maxCount:
            FrequentPatterns.append(item[0])

    return FrequentPatterns


def main():
    # read the data
    dat = np.loadtxt("dataset_9_8.txt", dtype=("str"))

    # Minimum Skew Problem
    # str_format = ""
    # for i in Skew(str(dat)):
    #     str_format = str_format + str(i) + " "
    # print(str_format)

    # Hamming Distance problem
    # print(HammingDistance(dat[0], dat[1]))

    # Approximate Pattern Matching Problem
    # str_format = ""

    # for i in ApproxPatternMatch(dat[1].strip('\n'), dat[0].strip('\n'), int(dat[2].strip('\n'))):
    #     str_format = str_format + str(i) + " "
    # print(str_format)

    # Count d problem

    # print(PatternCountD(dat[1].strip('\n'), dat[0].strip('\n'), int(dat[2].strip('\n'))))

    # print(Neighbours("CAA", 1))
    # ofile = open("outputfile.txt","w")

    # ofile.write('\n'.join(Neighbours("GACACCACCG", 3)))

    # ofile.close()

    # most frequent k-mer with up to d mismatches
    # str_format = ""

    # for i in FrequentWordsWithMismatch(dat[0].strip('\n'), int(dat[1].strip('\n')), int(dat[2].strip('\n'))):
    #     str_format = str_format + str(i) + " "
    # print(str_format)

    # Frequent Words with Mismatches and Reverse Complements Problem
    # str_format = ""

    # for i in FrequentWordsWithMismatchRC(dat[0].strip('\n'), int(dat[1].strip('\n')), int(dat[2].strip('\n'))):
    #     str_format = str_format + str(i) + " "
    # print(str_format)

if __name__ == "__main__":
    main()
