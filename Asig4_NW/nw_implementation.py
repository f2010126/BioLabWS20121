import fastaparser
import numpy as np
import pandas as pd

PAM_GAP = -8
BLOSUM62_GAP = -6


class NeedlemanWunsch:

    def __init__(self):
        super().__init__()

    def read_fasta_file(self, fasta_file):
        with open(fasta_file) as fasta_file:
            parser = fastaparser.Reader(fasta_file)
            for seq in parser:
                sequence = seq.sequence_as_string()

        return sequence

    def read_substitution_matrix(self, file_substitution_matrix):
        """
        Implement reading the scores file.
        It can be stored as a dictionary of example:
        scores[("A", "R")] = -1
        """
        submatrix = pd.read_csv(file_substitution_matrix,
                                sep=" ",
                                index_col=0,
                                header=0,
                                comment='#',
                                skipinitialspace=True)
        return submatrix

    def init_matrix(self, sequence1, sequence2, gap_cost):
        """
        Implement initialization of the matrix.
        Make sure you picked the right dimention and correctly initilized the first row and the first column.
        """
        n = len(sequence1) + 1
        m = len(sequence2) + 1
        matrix = np.zeros((n, m))
        # add values i * penalty
        for i in range(n):
            matrix[i][0] = i * gap_cost
        for j in range(m):
            matrix[0][j] = j * gap_cost
        return matrix

    def new_value_computation(self, char_seq1, char_seq2, gap_cost, substitution_matrix, diag_val, top_val, left_val):
        """
        Implement the computation of the value in the new cell.
        In this function we assume that we want to compute the value in the new cell in the matrix.
        Assume that the values "to the left", "to the top" and "top left" are already computed and provided
        as the input to the function. Also we know what characters in both sequences correspond to the given cell.
        """

        match = diag_val + substitution_matrix[char_seq1][char_seq2]
        cell_value = max(top_val + gap_cost, match, left_val + gap_cost)  # cell_value
        return cell_value

    def calculate_matrix(self, sequence1, sequence2, gap_opening_cost, substitution_cost):
        """
        Implement the step of complete computation of the matrix
        First initialize the matrix then fill it in from top to bottom.
        """
        matrix = self.init_matrix(seq1, seq2, gap_cost)
        for i in range(1, len(seq1)+1):

            for j in range(1, len(seq2)+1):
                top = matrix[i][j-1] # insert
                left = matrix[i-1][j] # delete
                diag = matrix[i-1][j-1]
                matrix[i][j] = self.new_value_computation(seq1[i-1],seq2[j-1],gap_cost,scoreMatrix,diag,top,left)

        return matrix


if __name__ == '__main__':

    algo = NeedlemanWunsch()
    # seq1 = algo.read_fasta_file("data/s1.fasta")
    # seq2 = algo.read_fasta_file("data/s2.fasta")
    gap_cost = BLOSUM62_GAP
    scoreMatrix = algo.read_substitution_matrix("data/blosum62.txt")
    seq1 = "AATC"
    seq2 = "AACGA"
    finished_matrix = algo.calculate_matrix(seq1, seq2, gap_cost, scoreMatrix)
    print(finished_matrix)
    #traceback = algo.traceback(finished_matrix,seq1,seq2,gap_cost,scoreMatrix)
    #print(traceback)
    # final_alignments = algo.alingment_build(traceback, seq1, seq2)
    # print(final_alignments)
    print("Done")
    # Tool to test o/p https://gtuckerkellogg.github.io/pairwise/demo/
