from math import log, sqrt
import numpy as np

from source.Parser import *

class GOR:
    def __init__(self, Parser):
        """
        GOR III implementation, divided in a training part and a prediction. The training creates the
        counters/frequency used afterwards to predict which structure an amino acid is part of according to the
        propensity of this amino acid R to be part of a structure S (I(S,R)) and its neighbourhood if it is part
        of this structure (I(S,R,Rj) with -8<=j<=8).
        :param Parser: DSSPParser object
        """
        self.freq_distance = np.zeros((3, 20, 20, 16))
        self.freq_struct = np.zeros(3)
        self.freq_acide = np.zeros((3, 20))
        for tup in Parser.get_trainset() :
            self.update_counters(tup[0], tup[1])


    def update_counters(self, sequence, structure):
        """
        Update the counters as asked : structure in which is the central amino acid, as well as its neighbourhood.
        :param sequence: str, amino acid sequence
        :param structure: str, sequence of the related structure
        :return: None, update the counters
        """
        for i in range(len(sequence)-1) :
            self.freq_struct[self.ret_struct(structure[i])] += 1
            self.freq_acide[self.ret_struct(structure[i])][self.ret_acid(sequence[i])] += 1
            for j in range(16):
                index = j - 8 + i
                if 0 <= index < len(sequence) and index != i:
                    self.freq_distance[self.ret_struct(structure[i])][self.ret_acid(sequence[index])][self.ret_acid(sequence[i])][j] += 1


    def predict(self, sequence):
        """
        Compute the probability of the three structures for each amino acid according to the counters and
        select the highest as the prediction (I(S,R) + somme(I(S,R,Rj) with -8<=j<=8 and j=!0))
        :param sequence: str, sequence to predict
        :return: str, sequence of the predicted structure
        """
        result = ""
        for i in range(len(sequence)):
            #result += self.GOR1(sequence[i])
            result += self.GOR3(sequence, i)
        return result


    def GOR1(self, acid):
        maxi = 0
        for stru in range(3) :
            res = self.p(acid, stru)
            if res > maxi :
                maxi = res
                index = stru
        return ['H', 'C', 'E'][index]


    def p(self, acid, i):
        other = (i + 1) % 3
        other_2 = (i + 2) % 3
        return log(self.freq_acide[i][self.ret_acid(acid)], 2) - log(self.freq_acide[other][self.ret_acid(acid)] + self.freq_acide[other_2][self.ret_acid(acid)], 2) \
              + log(self.freq_struct[other] + self.freq_struct[other_2], 2) - log(self.freq_struct[other], 2)


    def GOR3(self, sequence, i):
        maxi = -1
        for struct in range(3):
            res = self.p(sequence[i], struct) + sum(self.k(sequence, struct, i, m) for m in range(16))
            if res > maxi:
                maxi = res
                index = struct
        return ['H', 'C', 'E'][index]


    def k(self,sequence, struct, i, m):
        other = (struct + 1) % 3
        other_2 = (struct + 2) % 3
        acid = sequence[i]
        pos = i + m - 8
        if 0 <= pos < len(sequence) and m != 8:
            return log(self.freq_distance[struct][self.ret_acid(sequence[pos])][self.ret_acid(acid)][m], 2) - log(self.freq_distance[other][self.ret_acid(sequence[pos])][self.ret_acid(acid)][m] + \
                         self.freq_distance[other_2][self.ret_acid(sequence[pos])][self.ret_acid(acid)][m], 2) + log(self.freq_acide[other][self.ret_acid(acid)] + self.freq_acide[other_2][self.ret_acid(acid)],2) \
                   - log(self.freq_acide[struct][self.ret_acid(acid)],2)
        else :
            return 0

    def validate(self, test_set):
        """
        Specific function to test with a test set of variable length
        :param test_set: list of tuples of size 2 (str, str), respectively amino acid sequence and structure sequence
        :return: tuple of size 4 each of size 2((float, float), (float, float), (float, float), (float, float)),
        respectively mean Q3, MCC-H, MMC-E and MCC-C and their standard deviation, all rounded at 2 decimal places
        """
        Q3 = np.array([])
        MCC_H = np.array([])
        MCC_E = np.array([])
        MCC_C = np.array([])
        for tup in test_set:
            Q3 = np.append(Q3, self.Q3(tup[1], self.predict(tup[0])))
            H = self.MCC(tup[1], self.predict(tup[0]), 'H')
            E = self.MCC(tup[1], self.predict(tup[0]), 'E')
            C = self.MCC(tup[1], self.predict(tup[0]), 'C')
            if H is not None :  MCC_H = np.append(MCC_H, H)
            if E is not None:  MCC_E = np.append(MCC_E, E)
            if C is not None:  MCC_C = np.append(MCC_C, C)
        return (round(Q3.mean(),2), round(Q3.std(),2)), (round(MCC_H.mean(), 2), round(MCC_H.std(), 2)), (round(MCC_E.mean(), 2), round(MCC_E.std(), 2)), (round(MCC_C.mean(), 2), round(MCC_C.std(), 2))


    def Q3(self, real, predicted):
        """
        Q3 computation
        :param real: str, real structure sequence
        :param predicted: str, predicted structure sequence
        :return: float rounded at 2 decimal, percentage of good predictions overall
        """
        res = 0
        for i in range(len(real)) :
            if real[i] == predicted[i] : res += 1
        return round(res/len(real)*100,2)


    def MCC(self, real, predicted, structure):
        """
        Matthew coefficient correlation. Evaluates the performance of a classifier (here, we use the binary variation)
        Needs to have the True Positives, True Negatives, False Positives and False Negatives computed according to
        the structure given in parameter. Returns None in case of the denominator is non valid (division by 0)
        :param real: str, real structure sequence
        :param predicted: str, predicted structure sequence
        :param structure: str, structure you want to evaluate the MCC of
        :return: float rounded at 2 decimal, MCC value, or None if division error
        """
        TP = 0
        FP = 0
        FN = 0
        TN = 0
        for i in range(len(real)) :
            if real[i] == predicted[i] == structure:
                TP += 1
            if predicted[i] == structure and real[i] != structure :
                FP += 1
            if predicted[i] != structure and real[i] == structure :
                FN += 1
            if predicted[i] != structure and real[i] != structure :
                TN += 1
        root = (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
        if root != 0 :
            MCC = round(((TP * TN - FP * FN) / sqrt(root))*100,2)
        else :
            MCC = None
        return MCC


    def ret_acid(self, letter):
        acide = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        return acide.index(letter)


    def ret_struct(self, letter):
        struct = ['H', 'C', 'E']
        return struct.index(letter)


if __name__ == '__main__':
    fil = "D:/74ben/Documents/Programmation/PycharmProjects/project-3-secondary-structure-prediction-benjaminperraud/ressources/dataset/"

    parser = DSSPParser(fil)
    gor = GOR(parser)
    seq = 'RTDCYGNVNRIDTTGASCKTAK'
    real = 'CCCCCCCHHHCCCCCECHHHHH'
    predict = gor.predict(seq)