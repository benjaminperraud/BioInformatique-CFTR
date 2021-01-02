from math import *
import numpy as np


class Matrix:
    def __init__(self, nrows=0, ncols=0, value=None):
        self.matrix = np.full([nrows, ncols], value)
        self.row = nrows
        self.col = ncols

    def __getitem__(self, item):
        x, y = item
        return self.matrix[x][y]

    def get_num_cols(self):
        return len(self.matrix[0])

    def get_num_rows(self):
        return len(self.matrix)

    def get_max(self):
        """
        :return: tuple of size 3 of the row index, col index and value of the cell with the maximum value in the matrix
        If several occurences of the maximum value are found, the lowest indices are picked.
        """
        max = np.amax(self.matrix)
        index = np.unravel_index(self.matrix.argmax(), self.matrix.shape)
        return index[0], index[1], max

    def set_value(self, i, j, value):
        if 0 <= i < self.get_num_rows() and 0 <= j < self.get_num_cols():
            self.matrix[i][j] = value


class PSSM(Matrix):

    def __init__(self, multiple_alignment):
        super().__init__()
        self.PROB_AA = {'A': 0.0828, 'Q': 0.0394, 'L': 0.0967, 'S': 0.065,
                        'R': 0.0553, 'E': 0.0676, 'K': 0.0585, 'T': 0.0532, 'N': 0.0405,
                        'G': 0.0709, 'M': 0.0243, 'W': 0.0107, 'D': 0.0545, 'H': 0.0227,
                        'F': 0.0386, 'Y': 0.0291, 'C': 0.0136, 'I': 0.0599, 'P': 0.0468,
                        'V': 0.0687}
        self.create_mat(multiple_alignment)
        self.count_occurences(multiple_alignment)
        self.PWM()
        self.calculate_probabilities()

    def create_mat(self, file):
        file_submat = open(file, "r")
        self.acide = []
        self.col = 0
        self.nb_seq = 0
        for line in file_submat:
            if line[0] == '>':
                pass
            else:
                self.nb_seq += 1
                if self.col == 0: self.col = len(line) - 1
                for elem in line:
                    if elem not in self.acide and elem not in (' ', '-', '\n'):
                        self.acide.append(elem)
        self.row = len(self.acide)
        self.matrix = np.zeros([self.row, self.col])
        file_submat.close()

    def count_occurences(self, file):
        """
        Count the occurences of amino acids at each position
        """
        file_submat = open(file, "r")
        for line in file_submat:
            if line[0] == '>':
                pass
            else:
                i = 0
                for elem in line:
                    if elem in self.acide:
                        row = self.acide.index(elem)
                        self.set_value(row, i, self.matrix[row, i] + 1)
                    i += 1
        file_submat.close()

    def PWM(self):
        """
        Position-weight matrix is a matrix using the occurences of amino acids at specific position and
        normalize the counts according to the number of sequences to obtain the frequency of the amino acid at
        a position. We use here also pseudocount (beta) to avoid a null frequency at a position. The gaps are ignored.
        Corresponds to the q(u,a) formula
        """
        alpha = self.nb_seq - 1
        beta = sqrt(self.nb_seq)
        for i in range(self.row):
            for j in range(self.col):
                self.set_value(i, j, (alpha * self.matrix[i, j] / self.nb_seq + beta * self.PROB_AA[self.acide[i]]) / (alpha + beta))

    def calculate_probabilities(self):
        """
        To do the transition from PWM to PSSM, you have to transform the normalized frequency to a probability, with
        the help of the log-odds. This is the logarithm between the probability of an event and
        of this event occurring randomly. The random event here is the probability of finding the amino acid anywhere,
        the probability of the event is our observation of having this amino acid at this specific position.
        Corresponds to m(u,a) formula
        """
        for i in range(self.row):
            for j in range(self.col):
                self.set_value(i, j, round(log(self.matrix[i, j]/self.PROB_AA[self.acide[i]], 2), 2))

    def __getitem__(self, item):
        x, y = item
        try:
            return self.matrix[self.acide.index(x), y]
        except:
            return 0


if __name__ == '__main__':
    DIR = "D:/74ben/Documents/Programmation/PycharmProjects/Bio-Info2/ressources"
    pssm = PSSM(f"{DIR}/msaresults-MUSCLE.fasta")


