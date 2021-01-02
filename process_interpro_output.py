# To convert the interpro output to uniprot retrieve id list

FILE = "AAA_interpro.txt"
OUT = "AAA_uniprot.txt"

outfile = open(OUT,'w')

with open(FILE) as file:
    for line in file.readlines():
        line = line.strip().split('\t')
        id_acc,pos = line[0],line[-1]
        for p in pos.split(','):
            p = p.replace('..','-')
            outfile.write(f'{id_acc}[{p}]\n')