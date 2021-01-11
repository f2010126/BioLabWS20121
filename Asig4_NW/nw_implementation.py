import fastaparser # for reading the fasta file
import numpy as np
import pandas as pd


PAM_GAP = -8
BLOSUM62_GAP = -6
BLOSUM62_PATH = "data/blosum62.txt"
PAM_PATH = "data/pam250.txt"


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

        Args:
            file_substitution_matrix: file path with matrix scores

        Returns:
            dataframe with scores for easier lookup

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

        Args:
            sequence1: The first input sequence
            sequence2: The second i/p sequence
            gap_cost: gap openeing cost

        Returns:
            matrix with first row and column initialised with the gap penalty
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

        Args:
            char_seq1: character from seq 1
            char_seq2: character from seq 2
            gap_cost: linear gap cost
            substitution_matrix: scores to loop up
            diag_val: value from diagonal
            left_val: value from left
            top_val: value from top

        Returns:
            value of current cell


        """
        match = diag_val + substitution_matrix[char_seq1][char_seq2]
        cell_value = max(top_val + gap_cost, match, left_val + gap_cost)  # cell_value
        return cell_value

    def calculate_matrix(self, sequence1, sequence2, gap_opening_cost, substitution_cost):
        """
        Implement the step of complete computation of the matrix
        First initialize the matrix then fill it in from top to bottom.

        Args:
            sequence1: character from seq 1
            sequence2: character from seq 2
            gap_opening_cost: linear gap cost
            substitution_cost: scores to loop up

        Returns:
            completed matrix
        """
        matrix = self.init_matrix(sequence1, sequence2, gap_opening_cost)
        for i in range(1, len(sequence1)+1):
            for j in range(1, len(sequence2)+1):
                top = matrix[i][j-1]  # insert
                left = matrix[i-1][j]  # delete
                diag = matrix[i-1][j-1]
                matrix[i][j] = self.new_value_computation(sequence1[i-1], sequence2[j-1],
                                                          gap_opening_cost, substitution_cost,
                                                          diag, top, left)

        return matrix

    def traceback(self, matrix, sequence1, sequence2, gap_opening_cost, substitution_cost):
        """
        Implement the traceback part of the algorithm
        With the given matrix traceback which cells were taken to complete the path from
        the top left corner to the bottom right corner.

        Args:
            matrix: completed score matrix
            sequence1: character from seq 1
            sequence2: character from seq 2
            gap_opening_cost: linear gap cost
            substitution_cost: scores to loop up

        Returns:
            traceback to calculate optimal alignment
        """
        i = len(sequence1)
        j = len(sequence2)
        align_score = matrix[i][j]  # score in bottom right cell
        traceback = []
        while i > 0 and j > 0:
            score = matrix[i][j]
            diag = matrix[i-1][j-1]
            left = matrix[i][j-1]
            top = matrix[i-1][j]
            if score == (diag + substitution_cost[sequence1[i-1]][sequence2[j-1]]):
                traceback.append("diagonal")
                i -= 1
                j -= 1
            elif score == (top + gap_opening_cost):
                traceback.append("up")
                i -= 1
            elif score == (left + gap_opening_cost):
                traceback.append("left")
                j -= 1

        while i > 0:
            traceback.append("up")
            i -= 1

        while j > 0:
            traceback.append("left")
            j -= 1

        return traceback, align_score

    def alingment_build(self, traceback, sequence1, sequence2):
        """
        Implement the alingment creation.
        Given the traceback figure out which editing operations were used to create the alingment.

        Args:
            sequence1: character from seq 1
            sequence2: character from seq 2
            traceback: traceback to build alignment

        Returns:
            completed alignment strings
        """
        i = len(sequence1) - 1
        j = len(sequence2) - 1
        seq1_align = ""
        seq2_align = ""

        for step in traceback:
            if step == "diagonal":
                seq1_align += sequence1[i]
                seq2_align += sequence2[j]
                i -= 1
                j -= 1

            elif step == "up":
                seq1_align += sequence1[i]
                seq2_align += '-'
                i -= 1

            elif step == "left":
                seq1_align += '-'
                seq2_align += sequence2[j]
                j -= 1

        while j >= 0:
            seq1_align += '-'
            seq2_align += sequence2[j]
            j -= 1

        while i >= 0:
            seq2_align += '-'
            seq1_align += sequence1[i]
            i -= 1
        return seq1_align[::-1], seq2_align[::-1]


algo = NeedlemanWunsch()


def find_score(seq1, seq2, gap_cost, path):
    score_matrix = algo.read_substitution_matrix(path)
    finished_matrix = algo.calculate_matrix(seq1, seq2, gap_cost, score_matrix)
    traceback, alignment_score = algo.traceback(finished_matrix, seq1, seq2, gap_cost, score_matrix)
    align1, align2 = algo.alingment_build(traceback, seq1, seq2)
    print(f"Score: {alignment_score}")
    print(f"{align1}\n{align2}")


def test_small():
    find_score('ATC', 'GGAACT', BLOSUM62_GAP, BLOSUM62_PATH)


def test_s1s2():
    seq1 = algo.read_fasta_file("data/s1.fasta")
    seq2 = algo.read_fasta_file("data/s2.fasta")
    find_score(seq1, seq2, BLOSUM62_GAP, BLOSUM62_PATH)


def test_s1s3():
    seq1 = algo.read_fasta_file("data/s1.fasta")
    seq2 = algo.read_fasta_file("data/s3.fasta")
    find_score(seq1, seq2, BLOSUM62_GAP, BLOSUM62_PATH)

def test_s1s4():
    seq1 = algo.read_fasta_file("data/s1.fasta")
    seq2 = algo.read_fasta_file("data/s4.fasta")
    find_score(seq1, seq2, BLOSUM62_GAP, BLOSUM62_PATH)


if __name__ == '__main__':
    # test_small()  # test a small sequence
    test_s1s2()
    # Tool to test o/p https://gtuckerkellogg.github.io/pairwise/demo/
