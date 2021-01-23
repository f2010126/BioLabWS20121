import math  # need for infinity and Nan


def gotoh(fasta_file_1, fasta_file_2, cost_gap_open, file_substitution_matrix=None):
    """Put your code here"""
    algo = Gotoh(file_substitution_matrix)
    alignment_score, alignments = algo.run(fasta_file_1, fasta_file_2, cost_gap_open, file_substitution_matrix)
    return alignment_score, alignments


class Gotoh:

    def __init__(self, file_substitution_matrix):
        print("")
        self.score_matrix, self.lookup = self.set_up_substitution_matrix(file_substitution_matrix)

    def run(self, fasta_file_1, fasta_file_2, cost_gap_open, file_substitution_matrix=None):
        """Put your code here"""
        # seq1 = self.read_fasta_file(fasta_file_1)
        # seq2 = self.read_fasta_file(fasta_file_2)
        seq1 = 'CG'
        seq2 = 'CCGA'
        p_matrix = self.init_matrix_p(seq1, seq2)
        q_matrix = self.init_matrix_q(seq1, seq2)
        d_matrix = self.init_matrix_d(seq1, seq2, SCORES_DNA["alpha"], SCORES_DNA["beta"])

        return seq1, seq2  # alignment_score, alignments

    def read_fasta_file(self, fasta_file):
        """Implement reading the fasta file

        Args:
            fasta_file: file loctaion of sequence

        Returns:
              sequence
        """
        for line in open(fasta_file):
            li = line.strip()
            if not li.startswith(">"):
                return line.rstrip()  # sequence

    def set_up_substitution_matrix(self, file_substitution_matrix):
        """
        Args:
            file_substitution_matrix: location of lookup

        Returns:
            score matrix and lookup
        """
        scores = []
        for line in open(file_substitution_matrix):
            li = line.strip()
            if not li.startswith("#"):
                scores.append(line.rstrip().split())
        del scores[0]  # remove the first row of letters
        for row in scores:  # remove letter from each column
            del row[0]

        lookup = {"A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7, "H": 8,
                  "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, "S": 15, "T": 16,
                  "W": 17, "Y": 18, "V": 19}
        return scores, lookup

    def read_substitution_matrix(self, char1, char2):
        """
        Implement reading the scores file.
        It can be stored as a dictionary of example:
        scores[("A", "R")] = -1

        Args:
            char1: character from seq1
            char2: character from seq1

        Returns:
            score based on the lookup

        """
        if char1 in self.lookup and char2 in self.lookup:
            return self.score_matrix[self.lookup[char1]][self.lookup[char2]]
        else:
            return -8

    def init_matrix_d(self, seq_1, seq_2, cost_gap_open, cost_gap_extend):
        """
        Implement initialization of the matrix D
        Args:
            seq_1: first sequence
            seq_2: second sequence
            cost_gap_open:
            cost_gap_extend:
        """
        n = len(seq_1) + 1
        m = len(seq_2) + 1

        matrix_d = [[0 for i in range(m)] for j in range(n)]

        # add values open + i * extend
        for i in range(1, n):
            matrix_d[i][0] = self.costFunction(i, cost_gap_open, cost_gap_extend)
        for j in range(1, m):
            matrix_d[0][j] = self.costFunction(j, cost_gap_open, cost_gap_extend)
        return matrix_d

    def costFunction(self, i, open, extend):
        return open + i*extend

    def init_matrix_p(self, seq_1, seq_2):
        """
        Implement initialization of the matrix P
        Args:
            seq_1: first sequence
            seq_2: second sequence
        Returns:
            initialised P matrix
        """
        n = len(seq_1) + 1
        m = len(seq_2) + 1
        matrix_p = [[0 for i in range(m)] for j in range(n)]
        for i in range(1, n):
            matrix_p[i][0] = math.nan

        for j in range(1, m):
            matrix_p[0][j] = -math.inf
        matrix_p[0][0] = 'X'
        return matrix_p

    def init_matrix_q(self, seq_1, seq_2):
        """
        Implement initialization of the matrix Q
        Args:
            seq_1: first sequence
            seq_2: second sequence
        Returns:
            initialised P matrix
        """
        n = len(seq_1) + 1
        m = len(seq_2) + 1
        matrix_q = [[0 for i in range(m)] for j in range(n)]
        for i in range(1, n):
            matrix_q[i][0] = -math.inf

        for j in range(1, m):
            matrix_q[0][j] = math.nan
        matrix_q[0][0] = 'X'
        return matrix_q

    def complete_d_p_q_computation(self, seq_1, seq_2, cost_gap_open, cost_gap_extend, substitutions=None):
        """
        Implement the recursive computation of matrices D, P and Q
        """
        # TODO:


    """
    You are working with 3 matrices simultaneously.
    You can store your path as a list of cells.
    A cell can be a tuple: coordinates, matrix_name.
    And coordinates is a tuple of indexex i, j.

    Cell example: ((0, 2), "d")
    Path example: [((2, 4), 'd'), ((2, 4), 'q'), ((2, 3), 'q'), ((2, 2), 'd'), ((1, 1), 'd'), ((0, 0), 'd')]
    """

    def compute_all_tracebacks(self, seq1, seq2, d_matrix, p_matrix, q_matrix,
                               cost_gap_open, cost_gap_extend, substitution=None):
        """
        Implement a search for all possible paths from the bottom right corner to the top left.
        Implement 'find_all_previous' and check_complete first.

        """
        # TODO:
        return []  # all_paths

    def find_all_previous(self, cell, seq1, seq2, d_matrix, p_matrix, q_matrix,
                          cost_gap_open, cost_gap_extend, substitution=None):
        parent_cells = []
        """
        Implement a search for all possible previous cells.
        """
        # TODO:
        return parent_cells

    def check_complete(self, path):
        """
        Implement a function which checks if the traceback path is complete.
        """
        # TODO:
    def alignment(self, traceback_path, seq1, seq2):
        """
        Implement creation of the alignment with given traceback path and sequences1 and 2
        """
        # TODO:
        return '-', '-'  # alignment_seq1, alignment_seq2

    # HELPER FUNCTIONS FOR SELF CHECK

    def visualize_matrix(self, matrix):
        """
        Implement the visualization of a matrix.
        Can be used for self check
        """
        for line in matrix:
            print(line)

    def score_of_alignment(self, align_seq1, align_seq2, cost_gap_open,
                           cost_gap_extension, substitutions=None):
        """
        A nice helper function which computes the score of the given alignment.
        This is only used for the self check.
        Input example:
        --CG
        AACA
        """
        # TODO:
        return 0  # score


SCORES_DNA = {'match': 1,
              'mismatch': -1,
              'alpha': -3,
              'beta': -1}

SCORES_PRO = {'alpha': -11,
              'beta': -1}
if __name__ == '__main__':

    fasta1 = "data/s1.fasta"
    fasta2 = "data/s2.fasta"
    pam_file = "data/pam250.txt"
    blosum_file = "data/blosum62.txt"
    gap_open = SCORES_DNA
    gotoh(fasta1, fasta2, SCORES_DNA, pam_file)
    # Tool to test o/p http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh
