import numpy as np

class Matrix:
    def __init__(self, nrows = 0, ncols = 0, value = None):
        self.matrix = np.full( [nrows, ncols], value)
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


class SubstitutionMatrix(Matrix):
    """
    Format : free
    """
    def __init__(self, file):
        super().__init__()
        self.parse_file(file)

    def parse_file(self, file):
        file_submat = open(file, "r")
        first = True
        row = 0
        col = 0
        for line in file_submat :
            if not first :
                liste = line.split()[1:-4]                  # on exclut B, Z, X et *
                for l in liste :
                    self.set_value(row, col, int(l))
                    col += 1
                col = 0
                row += 1
            if line[0] != "#" and first :
                first = False
                self.abc = line.split()
                self.matrix = np.zeros((len(self.abc) -4 , len(self.abc) -4))  # on exclut B, Z, X et *
        file_submat.close()

    def __getitem__(self, item):
        x, y = item
        try:
            return self.matrix[self.abc.index(x), self.abc.index(y)]
        except:
            return 0


## Proposition of implementation, but free for you to use it or not, not included in tests
class ScoreMatrix(Matrix):
    """
    Format : list of list of float values
    """
    def __init__(self, seq1, seq2, I, E, local=False, type = "S"):
        super().__init__()


if __name__ == '__main__':
    submat = SubstitutionMatrix("E:/74ben/Documents/Programmation/PycharmProjects/BioInfo-ADN/ressources/pam120.txt")
    mat = Matrix(4,4)