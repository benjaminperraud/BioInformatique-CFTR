class DSSPParser:

    def __init__(self, directory, trainset_size=3000):
        """
        Read the DSSP files and retrieve the amino acid sequences of the chains and their structures
        Create a list of the sequences and structures for the training set and a list for the testing set
        :param directory: str, destination folder where we can find the dssp files
        :param trainset_size: number of proteins in the training set vs the proteins in the testing set
        """
        file = "ressources/dataset/CATH_info.txt"
        self.seq = []
        #cath = open(directory + "CATH_info.txt", "r")
        cath = open(file, "r")
        for line in cath:
            #self.parse(directory + "dssp/" + line.split()[0])
            self.parse(directory + '/' + line.split()[0])

        self.train = self.seq[:trainset_size]
        self.test = self.seq[trainset_size:]
        cath.close()

    def parse(self, file):
        """
        Extract the information for each identifiers in the CATH info file. Each identifier contains the id of the
        protein, as well as the chain of the desired sequence. The amino acid X, Z, B will be ignored, and all
        lower caps amino acids are converted to Cystein (C). The structures are converted into C (C, T, S and " "),
        H (H, G, I) or E (E and B).
        :return: None, but create a list of tuples of size 2 (str, str), respectively the amino acid sequence and the
        sequence of the structure
        """
        self.acid = ""
        self.struct = ""
        dssp = open(file[:-1]+".dssp", "r")
        group = file[-1]
        i = 0
        for line in dssp :
            if i >= 28 :
                l = line.split()
                if line[11] == group :
                    if all(chr.isdigit() for chr in l[1]) :
                        if line[13] not in ('*', '!', 'B', 'X', 'Z') :
                            if not line[13].isupper():
                                char = 'C'
                            else :
                                char = line[13]
                            self.acid += char
                            if line[16] in ('C', 'T', 'S', " "):
                                self.struct += 'C'
                            if line[16] in ('H', 'G', 'I'):
                                self.struct += 'H'
                            if line[16] in ('E', 'B'):
                                self.struct += 'E'
            i += 1
        self.seq.append((self.acid, self.struct))
        dssp.close()

    def get_trainset(self):
        return self.train

    def get_testset(self):
        return self.test


if __name__ == '__main__':
    fil = "D:/74ben/Documents/Programmation/PycharmProjects/project-3-secondary-structure-prediction-benjaminperraud/ressources/dataset/"
    parse = DSSPParser(fil)
    print(parse.seq[115])
