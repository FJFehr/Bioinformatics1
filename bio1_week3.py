

from bio1_week1 import NumberToPattern
from bio1_week2 import Neighbours
from bio1_week2 import HammingDistance
import numpy as np


def MotifEnumeration(Dna, k, d):

    """
    Dna is an array of strings
    """

    Patterns = []
    kmer_array = []
    neighbour_array = []
    temp_array1 = []
    temp_array2 = []

    # split up DNA strings into kmers
    for i in range(len(Dna)):
        for j in range(len(Dna[i])-k+1):

            # Populate kmer temp array
            temp_array1.insert(j, Dna[i][j:j+k])

            # For the kmer above do the neighbours
            neighbours_temp = Neighbours(Dna[i][j:j+k], d)

            for n in neighbours_temp:
                temp_array2.append(n)

        # Populate dna nested array
        kmer_array.insert(i, temp_array1)

        neighbour_array.insert(i, temp_array2)

        temp_array1 = []
        temp_array2 = []

    Patterns = set(neighbour_array[0]).intersection(*neighbour_array)

    return Patterns


def DistanceBetweenPatternAndStrings(Pattern, Dna):

    '''
    Dna is an array of Dna strings
    Returns the distance between a Patterns and a List of Dnas
    Needed for MedianString

    '''

    k = len(Pattern)
    Distance = 0

    for Text in Dna:
        HamDist = k  # Very large number

        for i in range(len(Text) - k + 1):

            currentPattern = Text[i: i + k]

            if HamDist > HammingDistance(currentPattern, Pattern):
                HamDist = HammingDistance(currentPattern, Pattern)

        Distance = Distance + HamDist

    return Distance


def MedianString(Dna, k):

    """
    Returns the median string of an array of DNA sequences
    """
    distance = 10000

    for i in range(4**k-1):

        Pattern = NumberToPattern(i, k)

        if distance > DistanceBetweenPatternAndStrings(Pattern, Dna):
            distance = DistanceBetweenPatternAndStrings(Pattern, Dna)
            median = Pattern

    return median


def ProfileMostProbable(Dna, k, MatrixProfile):
    """
    Dna is a string
    k is an int
    MatrixProifile is a numpy matrix
    """

    # first get a list of all k_mers
    kmerList = list()
    for i in range(len(Dna)-k +1):
        kmerList.append(Dna[i:i+k])

    # Calc probs of all kmers in list above
    kmerProbDict = dict()

    for kmer in kmerList:

        # prob of current kmer
        currentProb = 1

        for i in range(len(kmer)):

            if kmer[i] == "A":  # row 1
                currentProb = currentProb * MatrixProfile[0, i]
            if kmer[i] == "C":  # row 2
                currentProb = currentProb * MatrixProfile[1, i]
            if kmer[i] == "G":  # row 3
                currentProb = currentProb * MatrixProfile[2, i]
            if kmer[i] == "T":  # row 4
                currentProb = currentProb * MatrixProfile[3, i]

        kmerProbDict[kmer] = currentProb

    # Get the max value
    return max(kmerProbDict, key=kmerProbDict.get)


def Profile(Motifs, laplace_bool):
    """
    Motifs are an array of strings length txk
    Meaning t strings of length k
    Returns a Numpy matrix of size txk
    """

    # we assume all motifs are of the same length k
    k = len(Motifs[0])
    # Length of all motifs
    t = len(Motifs)

    counts_dict = {"A": 0, "C": 0, "G": 0, "T": 0}

    # Create the empty numpy matrix
    current_profile = np.zeros((4, k))

    # Add 1 to all entries if laplace is wanted
    if laplace_bool:
        current_profile += 1

    # Make profiles
    for i in range(k):
        for j in range(t):
            # Get whole column in one list
            counts_dict[Motifs[j][i]] = counts_dict[Motifs[j][i]] + 1

        # take into account laplace
        if laplace_bool:
            current_max = t + 4
        else:
            current_max = t

        # Place counts or profiles into current place
        counter = 0
        for val in counts_dict:
            current_profile[counter, i] = (counts_dict[val] + current_profile[counter, i])/current_max
            counter += 1

        # reset the dictionary
        counts_dict = {"A": 0, "C": 0, "G": 0, "T": 0}

    return current_profile


def Entropy(Motifs):
    """
    Motifs are an array of strings length txk
    Meaning t strings of length k
    returns a scaler score value
    """

    return Motifs


def EvaluateKmer(Profile, SingleDna):
    """
    PR(Pattern = SingleDNA | Profile)
    Evaluate the probability of a single kmer
    Profile is a Numpy matrix of size txk
    SingleDna is a string of length 1xk (Kmer)
    """

    index_dict = {"A": 0, "C": 1, "G": 2, "T": 3}

    current_prob = 1
    counter = 0

    for i in SingleDna:

        row_index = index_dict[i]

        current_prob = current_prob * Profile[row_index, counter]

        counter += 1

    return current_prob


def Motifs(Profile, Dna):
    """
    Collection of most probable k-mers formed by the Profile most probable
    Profile is a Numpy matrix of size txk
    DNA is an array of strings length txn
    t is number of strands and n is the length.
    returns an array of motif strings
    """

    k = Profile.shape[1]
    t = Profile.shape[0]

    motifs_array = []

    # iterate through DNA strands
    for strand in Dna:

        kmer_prob_dict = {}

        for i in range(len(strand) - k + 1):
            # get current kmer
            current_kmer = strand[i:i+k]
            # evaluate the prob
            current_kmer_prob = EvaluateKmer(Profile, current_kmer)
            # append it to dictionary
            kmer_prob_dict[current_kmer] = current_kmer_prob

        #get key with max value and append to motifs array
        motifs_array.append(max(kmer_prob_dict, key=kmer_prob_dict.get))

    return motifs_array


def Score(Motifs):
    """
    Motifs are an array of strings length txk
    Meaning t strings of length k
    returns a scaler score value
    """
    # we assume all motifs are of the same length k
    k = len(Motifs[0])
    # Length of all motifs
    t = len(Motifs)

    # Holds the first column
    counts_dict = {"A": 0, "C": 0, "G": 0, "T": 0}
    current_score = 0

    for i in range(k):
        for j in range(t):
            # Get whole column in one list
            counts_dict[Motifs[j][i]] = counts_dict[Motifs[j][i]] + 1

        # find max in counts dict at each iteration and then go t-max
        current_score += t - max(counts_dict.values())

        # reset the dictionary
        counts_dict = {"A": 0, "C": 0, "G": 0, "T": 0}

    return current_score


def GreedyMotifSearch(Dna, k, t):

    bestMotifs = list()

    # Randomly assign the first k of each DNA to be best
    for strand in Dna:
        bestMotifs.append(strand[0:k])

    baseStrand = Dna[0]
    otherStrands = Dna[1:]

    for i in range(len(baseStrand) - k + 1):
        currentMotif = baseStrand[i:i+k]

        for strand in otherStrands:
            #We need to create a profile matrix now...
            profileMatrix = Profile(motifs)




    return baseStrand, otherStrands, bestMotifs


def main():
    # read the data
    # dat = np.loadtxt("test_dataset1.txt", dtype=("str"))

    # with open('dataset_156_8.txt', 'r') as myfile:
    #     dat = myfile.read().replace('\n', ' ')
    # myfile.close()

    # dat = dat.split()

    # new =[]
    # for i in dat[2:]:
    #     new.append(i.strip('\n'))

    # str_format = ""
    # res = sorted(MotifEnumeration(new, int(dat[0].strip('\n')), int(dat[1].strip('\n'))))
    # for i in res:
    #     str_format = str_format + str(i) + " "
    # print(str_format)

    # Test DistanceBetweenPatternAndStrings.

    # with open('dataset_5164_1.txt', 'r') as myfile:
    #     dat = myfile.read().replace('\n', ' ')
    # myfile.close()

    # dat = dat.split()

    # print(DistanceBetweenPatternAndStrings(dat[0], dat[1:]))

    # Median String
    # with open('dataset_158_9.txt', 'r') as myfile:
    #     dat = myfile.read().replace('\n', ' ')
    # myfile.close()

    # dat = dat.split()
    # print(MedianString(dat[1:], int(dat[0])))

    # Profile-most Probable k-mer Problem

    # with open('dataset_159_3.txt', 'r') as myfile:
    #     dat = myfile.read().replace('\n', ' ')
    # myfile.close()

    # dat = dat.split()

    # Dna = dat[0]

    # k = int(dat[1])

    # # Getting it to be soek matix
    # profile = np.mat(dat[2:])
    # profile = profile.reshape(4, k)
    # profile = profile.astype(np.float)

    # print(ProfileMostProbable(Dna, k, profile))

    # Implement GreedyMotifSearch

    # with open('test_dataset.txt', 'r') as myfile:
    #     dat = myfile.read().replace('\n', ' ')
    # myfile.close()

    # dat = dat.split()

    # print(GreedyMotifSearch(dat[2:], int(dat[0]), int(dat[1])))
    # Motifs = [
    #     "TCGGGGGTTTTT",
    #     "CCGGTGACTTAC",
    #     "ACGGGGATTTTC",
    #     "TTGGGGACTTTT",
    #     "AAGGGGACTTCC",
    #     "TTGGGGACTTCC",
    #     "TCGGGGATTCAT",
    #     "TCGGGGATTCCT",
    #     "TAGGGGAACTAC",
    #     "TCGGGTATAACC"
    # ]
    # print(Profile(Motifs2, True))



    current_profile = np.matrix([[0.8, 0, 0, 0.2], [0, 0.6, 0.2, 0], [0.2, 0.2, 0.8, 0], [0, 0.2, 0, 0.8]])
    current_dna = ["TTACCTTAAC", "GATGTCTGTC", "ACGGCGTTAG", "CCCTAACGAG", "CGTCAGAGGT"]

    #print(EvaluateKmer(current_profile, "ACCT"))

    print(Motifs(current_profile, current_dna))

if __name__ == "__main__":
    main()
