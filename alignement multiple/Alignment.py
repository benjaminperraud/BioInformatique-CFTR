from source.Matrix import *
from abc import ABC, abstractmethod


class Alignment(ABC):
    def __init__(self, seq1, seq2, I, E, submat, k):
        self.seq1 = seq1
        self.seq2 = seq2
        self.I = I
        self.E = E
        self.submat = submat
        self.mat_S = Matrix(len(seq1) + 1, len(seq2) + 1)  # matrice de score
        self.mat_V = Matrix(len(seq1) + 1, len(seq2) + 1)  # matrice de gap horizontale
        self.mat_W = Matrix(len(seq1) + 1, len(seq2) + 1)  # matrice de gap verticale
        self.k = k
        self.solution = []  # list of all solution for alignement

    @abstractmethod
    def backtrack(self, i, j, ali1, ali2):
        """
        Use the Score Matrices to find the path taken and form the two strings of the alignment
        :param ali1: alignement 1
        :param ali2: alignement 2
        :param i: int, row index of current position
        :param j: int, col index of current position
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """

    @abstractmethod
    def run(self):
        """
        Run the alignment algorithm according to the parameters
        :return:
        """
        pass

    def compute_scores(self, alignment):
        """
        Compute the identity, similarity and number of gaps of an alignment
        :param alignment: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2
        with the inserted gaps and score of the alignment and the positions of beginning and end of the alignments sequences
        :return: list of tuples of 3 floats (rounded two decimal places) respectively identity, similarity and gaps rates (in %)
        """
        sol = []
        for ali in alignment:
            seq1 = ali[0]
            seq2 = ali[1]
            identity = 0  # identity score
            similarity = 0  # similarity score
            gaps = 0  # gaps rate
            for i in range(len(seq1)):
                if seq1[i] == seq2[i]:
                    identity += 1
                if seq1[i] == "-":
                    gaps += 1
                if seq2[i] == "-":
                    gaps += 1
                if self.submat[seq1[i], seq2[i]] > 0:
                    similarity += 1
            sol.append((round(100 * identity / len(seq1), 2), round(100 * similarity / len(seq1), 2),
                        round(gaps / (0.02 * len(seq1)), 2)))
        return sol

    def get_solution(self):
        """
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """
        return self.solution


class SmithWaterman(Alignment):
    def __init__(self, seq1, seq2, I, E, submat, k):
        """
        Local alignment algorithm
        :param seq1: str, first sequence of amino acids to align
        :param seq2: str, second sequence of amino acids to align
        :param I: float, initial gap value
        :param E: float, extension gap value
        :param submat: SubstitutionMatrix object
        :param k: int, maximum number of solutions to find
        """
        super().__init__(seq1, seq2, I, E, submat, k)
        self.avoid = Matrix(len(seq1)+1, len(seq2)+1, 1).matrix   # garde en mémoire les indices qui sont mis à 0 (1 si pas encore visité)


    def backtrack(self, i, j, ali1, ali2):
        """
        Use the Score Matrices to find the path taken and form the two strings of the alignment
        :param ali1: alignement 1
        :param ali2: alignement 2
        :param i: int, row index of current position
        :param j: int, col index of current position
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """
        if self.mat_S.matrix[i,j] == 0:
            ali1.reverse()
            ali2.reverse()
            ali1 = "".join(ali1)
            ali2 = "".join(ali2)
            return ali1, ali2, self.start[2], (i+1,j+1), self.start[:-1]
        if self.mat_S[i, j] == self.mat_S[i - 1, j - 1] + self.submat[self.seq1[i - 1],self.seq2[j - 1]]: # match
            ali2.append(self.seq2[j - 1])
            ali1.append(self.seq1[i - 1])
            self.zero(i,j)
            return self.backtrack(i - 1, j - 1, ali1, ali2)
        if self.mat_S[i, j] == self.mat_V[i, j] :  # gap séquence horizontale
            ali2.append('-')
            ali1.append(self.seq1[i - 1])
            self.zero(i, j)
            return self.backtrack(i - 1, j, ali1, ali2)
        if self.mat_S[i, j] == self.mat_W[i, j] :  # gap séquence verticale
            ali1.append('-')
            ali2.append(self.seq2[j - 1])
            self.zero(i, j)
            return self.backtrack(i, j - 1, ali1, ali2)


    def zero(self, i, j):
        """
        Put all matrices at index i, j at 0
        :param i: int, row index of current position
        :param j: int, col index of current position
        :return:
        """
        self.mat_S.matrix[i, j] = 0
        self.mat_V.matrix[i, j] = 0
        self.mat_W.matrix[i, j] = 0
        self.avoid[i, j] = 0


    def run(self):
        """
        Run the alignment algorithm according to the parameters
        :return:
        """
        self.init_mat()
        for i in range(self.k):
            ali1 = []
            ali2 = []
            self.start = self.mat_S.get_max()
            self.solution.append(self.backtrack(self.start[0], self.start[1], ali1, ali2))
            self.recalculate(self.solution[i][3][0], self.solution[i][3][1])     # recalcul de la matrice depuis les indices d'arrivées du dernier alignement


    def init_mat(self):
        """
        Initialize matrices S, V and W
        :return:
        """
        # matrice Score
        self.mat_S.matrix.fill(0)
        # matrice V
        self.mat_V.matrix.fill(0)
        for i in range(0, self.mat_V.col):
            self.mat_V.set_value(0, i, -np.inf)
        # matrice W
        self.mat_W.matrix.fill(0)
        for j in range(0, self.mat_W.row):
            self.mat_W.set_value(j, 0, -np.inf)
        # Remplit les 3 matrices S, V et W
        self.recalculate(1,1)


    def recalculate(self, x, y):
        """
        The path taken must be erased (values put to 0 and all the values in the matrix below the last cell of the path
        must be recomputed. The values at 0 must stay at 0.
        :param x: row index of recalcul start
        :param y: col index of recalcul start
        :return: None, but the ScoreMatrix S has been modified accordingly
        """
        for i in range(x, self.mat_S.row):
            for j in range(y, self.mat_S.col):
                if self.avoid[i,j] != 0:
                    v = max(self.mat_S[i - 1, j] - self.I, self.mat_V[i - 1, j] - self.E)
                    self.mat_V.set_value(i, j, v if v > 0 else 0)
                    w = max(self.mat_S[i, j - 1] - self.I, self.mat_W[i, j - 1] - self.E)
                    self.mat_W.set_value(i, j, w if w > 0 else 0)
                    s = max(self.mat_S[i - 1, j - 1] + self.submat[self.seq1[i - 1], self.seq2[j - 1]], self.mat_W[i, j], self.mat_V[i, j])
                    self.mat_S.set_value(i, j, s if s > 0 else 0)


class SmithWatermanProfile(SmithWaterman):
    def __init__(self, seq, I, E, PSSM, k):
        """
        Local alignment algorithm
        :param seq: str, sequence of amino acids to find matches for the profile matrix
        :param I: float, initial gap value
        :param E: float, extension gap value
        :param PSSM: PSSM object
        :param k: int, maximum number of solutions to find
        """
        super().__init__(seq, ' ' * (PSSM.col), I, E, PSSM, k)  # bonne longueur de la seq2 pour la construction des matrices

    def run(self):
        """
        Run the alignment algorithm according to the parameters
        :return:
        """
        self.init_mat()
        for i in range(self.k):
            ali1 = []
            ali2 = []
            start = self.mat_S.get_max()
            back = self.backtrack(start[0], start[1], ali1, ali2)
            self.solution.append((back[0], round(start[2], 2), (back[1], start[0])))
            self.recalculate(back[1], back[2])  # recalcul de la matrice depuis les indices d'arrivées du dernier alignement


    def backtrack(self, i, j, ali1, ali2):
        """
        Use the Score Matrices to find the path taken and form the two strings of the alignment
        :param ali1: alignement 1
        :param ali2: alignement 2
        :param i: int, row index of current position
        :param j: int, col index of current position
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """
        if self.mat_S.matrix[i, j] == 0:
            ali1.reverse()
            ali1 = "".join(ali1)
            return ali1, i + 1, j + 1
        if self.mat_S[i, j] == self.mat_S[i - 1, j - 1] + self.submat[self.seq1[i-1], j-1]:  # match
            ali1.append(self.seq1[i - 1])
            self.zero(i, j)
            return self.backtrack(i - 1, j - 1, ali1, ali2)
        if self.mat_S[i, j] == self.mat_V[i, j]:  # gap séquence horizontale
            ali1.append(self.seq1[i - 1])
            self.zero(i, j)
            return self.backtrack(i - 1, j, ali1, ali2)
        if self.mat_S[i, j] == self.mat_W[i, j]:  # gap séquence verticale
            ali1.append('-')
            self.zero(i, j)
            return self.backtrack(i, j - 1, ali1, ali2)

    def recalculate(self, x, y):
        """
        The path taken must be erased (values put to 0 and all the values in the matrix below the last cell of the path
        must be recomputed. The values at 0 must stay at 0.
        :param x: row index of recalcul start
        :param y: col index of recalcul start
        :return: None, but the ScoreMatrix S has been modified accordingly
        """
        for i in range(x, self.mat_S.row):
            for j in range(y, self.mat_S.col):
                if self.avoid[i, j] != 0:
                    v = max(self.mat_S[i - 1, j] - self.I, self.mat_V[i - 1, j] - self.E)
                    self.mat_V.set_value(i, j, v if v > 0 else 0)
                    w = max(self.mat_S[i, j - 1] - self.I, self.mat_W[i, j - 1] - self.E)
                    self.mat_W.set_value(i, j, w if w > 0 else 0)
                    s = max(self.mat_S[i - 1, j - 1] + self.submat[self.seq1[i-1], j-1],
                            self.mat_W[i, j], self.mat_V[i, j])
                    self.mat_S.set_value(i, j, s if s > 0 else 0)


if __name__ == '__main__':
    DIR = "D:/74ben/Documents/Programmation/PycharmProjects/Bio-Info2/ressources"
    pssm = PSSM(f"{DIR}/msaresults-MUSCLE.fasta")

    seq = 'MADEEKLPPGWEKRMSRSSGRVYYFNHITNASQWERPSGNSSSGGKNGQGEPARVRCSHLLVKHSQSRRPSSWRQEKITRTKEEALELINGYIQKIKSGEEDFESLASQFSDCSSAKARGDLGAFSRGQMQKPFEDASFALRTGEMSGPVFTDSGIHIILRTE'
    align = SmithWatermanProfile(seq, 1, 1, pssm, 2)  # only linear test
    align.run()
    print(align.get_solution())