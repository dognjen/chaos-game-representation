import numpy
import matplotlib.cm as cm
import pylab

k=4             # number of combinations
s=2**k          # size of the matrix

bases = ["A", "C", "G", "T"]

def cartesian(nucleotides):
    """ Create a cartesian product from variable nucleotides. """
    if not nucleotides:
        yield ""
    else:
        for a in nucleotides[0]:
            for b in cartesian(nucleotides[1:]):
                yield a + b


def count(sequence, sub):
    """ Method counts substrings in a given sequence. Moving through the sequence by one nucleotide."""
    counter = 0
    for i in range(len(sequence) - len(sub) + 1):
        if sequence[i:i + len(sub)] == sub:
            counter += 1
    return counter


def k_mer_probability(sequence, bases):
    """ Method calculates the probability of k-mer in a given sequence. """
    probabilities = {}
    for b in bases:
        probabilities[b] = count(sequence, b) / float(len(sequence) - len(b) + 1)
    return probabilities


def fill_matrix_with_probabilities(bases, depth, row, column, nuk):
    """ This recursion fills the matrix of nucleotides with probabilities for chaos game representation.
     Nucleotide A represents upper left corner, nucleotide C represents bottom left corner, G represents bottom right,
     while T represents upper right corner of the image. """

    if depth == 0:
        mat[column[0]][row[0]] = (ver[nuk])
        return []
    else:
        for b in bases:
            if b in ['A', 'C']:
                r = row[:len(row)//2]
            else:
                r = row[len(row)//2:]
            if b in ['A', 'T']:
                c = column[:len(column)//2]
            else:
                c = column[len(column)//2:]
            fill_matrix_with_probabilities(bases, depth - 1, r, c, nuk + b)
        return []


if __name__ == "__main__":
    fasta = open("./NC_012920.1.fasta")
    fasta = fasta.read()
    sequence = "".join(fasta.split("\n")[1:])

    # remove unknown nucleotides
    sequence = sequence.replace("N", "")

    combinations = list(cartesian([bases] * k))
    ver = k_mer_probability(sequence, combinations)

    print(ver,"\n")
    print(fasta,"\n")

    mat = numpy.array([[0.0 for i in range(s)] for j in range(s)], float)

    for i in range(len(combinations)):
        print(i,i//s,i%s,ver[combinations[i]]*10)
        mat[i//s, i%s] = ver[combinations[i]]*10

    fill_matrix_with_probabilities(bases, k, [int(i) for i in range(s)], [int(i) for i in range(s)], '')

    pylab.imshow(mat, extent=[0, s, 0, s], interpolation='nearest', cmap=cm.gray)
    novix = [a + 0.4 for a in range(s)]
    noviy = [a + 0.4 for a in range(s)]
    imena_x = [str(i) for i in range(s)]
    imena_y = [str(i) for i in range(s-1, -1, -1)]
    pylab.xticks(novix, imena_x)
    pylab.yticks(noviy, imena_y)
    pylab.show()




