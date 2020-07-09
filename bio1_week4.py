

from bio1_week3 import Profile
from bio1_week3 import Score
from bio1_week3 import Motifs
from random import randint
import numpy as np


def RandomizedMotifSearch(Dna, k, t):

    """
    Dna is an array of DNA strings
    k is the length of the kmer
    t is the number of motifs (and rows in Dna)
    """

    # Starting choice

    actual_best_motifs = []
    for strand in Dna:
        # generate random number
        rand = randint(0, len(strand) - k)

        # append random motif to starting array
        actual_best_motifs.append(strand[rand:rand+k])

    # start the iterations
    for i in range(100):

        motif_array = []
        # Choose random starting motifs
        for strand in Dna:
            # generate random number
            rand = randint(0, len(strand) - k)
            # append random motif to starting array
            motif_array.append(strand[rand:rand+k])

        current_best_motifs = motif_array

        while True:

            current_profile = Profile(current_best_motifs, True)
            current_motifs = Motifs(current_profile, Dna)

            if Score(current_motifs) < Score(current_best_motifs):
                current_best_motifs = current_motifs

            else:
                break

        if Score(current_best_motifs) < Score(actual_best_motifs):
            actual_best_motifs = current_best_motifs

    return actual_best_motifs

def ProfileRandomlyGen(Dna, k, MatrixProfile):
    """
    Generates a single weighted profile from
    Dna is a string
    k is an int
    MatrixProifile is a numpy matrix
    """
    # We loop through dnastrand to make kmers and find PR(Kmer|profileLessOne)

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

    # Now we scale it to make a proper PMF
    sum_p = sum(kmerProbDict.values())
    for i in kmerProbDict:
        kmerProbDict[i] = float(kmerProbDict[i]/sum_p)

    # Proper PMF
    return kmerProbDict

def GibbsSampler(Dna, k, t, N):

    """
    Dna is an array of DNA strings
    k is the length of the kmer
    t is the number of motifs (and rows in Dna)
    N is iterations
    """

    # Starting choice

    motifs = RandomizedMotifSearch(Dna, k, t)

    # motifs = []
    # for strand in Dna:
    #     # generate random number
    #     rand = randint(0, len(strand) - k)

    #     # append random motif to starting array
    #     motifs.append(strand[rand:rand+k])

    bestMotifs = motifs

    # Start iterations

    for j in range(N):

        # rand number for random leave out
        rand = randint(0,t-1)

        motifsLessOne = [x for i, x in enumerate(motifs) if i!=rand]

        # Get profile matrix of motifsLessOne

        profileLessOne = Profile(motifsLessOne, True)

        # we ultimately need to replace the motif[i] with something from dnaStrand[i]

        # This is our PMF
        missingProfile = ProfileRandomlyGen(Dna[rand], k, profileLessOne)

        # Toss the coin and get a motif
        new_motif = np.random.choice(a=list(missingProfile.keys()), p=list(missingProfile.values()))

        # Add the new motif in
        motifsLessOne.insert(rand, new_motif)

        motifs = motifsLessOne

        if Score(motifs) < Score(bestMotifs):
            bestMotifs = motifs

    return bestMotifs


def main():
    with open('test_dataset.txt', 'r') as myfile:
        dat = myfile.read().replace('\n', ' ')
    myfile.close()

    dat = dat.split()

    # # RANDOMISED MOTIF SEARCH
    # array = RandomizedMotifSearch(dat[2:], int(dat[0]), int(dat[1]))

    # for i in array:
    #     print(i)

    # GIBBS SAMPLER

    array = GibbsSampler(dat[3:], int(dat[0]), int(dat[1]), int(dat[2]))

    for i in array:
        print(i)


if __name__ == "__main__":
    main()