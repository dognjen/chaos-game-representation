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

	# Create all possible k combinations of base nucleotides
    combinations = list(cartesian([bases] * k))
	
	# Calculate the probability (k-mer) of all combinations in a given sequence
    ver = k_mer_probability(sequence, combinations)

    print(ver,"\n")
    print(fasta,"\n")

	# create a matrix in which to represent the combinations
    mat = numpy.array([[0.0 for i in range(s)] for j in range(s)], float)

	# fill the matrix with probabilities, as you can tell by the function name
    fill_matrix_with_probabilities(bases, k, [int(i) for i in range(s)], [int(i) for i in range(s)], '')

	# visualization, 
    pylab.imshow(mat, extent=[0, s, 0, s], interpolation='nearest', cmap=cm.gray)
    x_os = [a + 0.4 for a in range(s)]
    y_os = [a + 0.4 for a in range(s)]
    label_x = [str(i) for i in range(s)]
    label_y = [str(i) for i in range(s-1, -1, -1)]
    pylab.xticks(x_os, label_x)
    pylab.yticks(y_os, label_y)
    pylab.show()




