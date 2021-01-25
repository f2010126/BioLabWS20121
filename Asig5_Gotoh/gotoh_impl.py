import math  # need for infinity and Nan


def gotoh(fasta_file_1, fasta_file_2, scores, is_dna=False, file_substitution_matrix=None):
    """Put your code here"""
    algo = Gotoh(file_substitution_matrix, scores, is_dna)
    alignment_score, alignments = algo.run(fasta_file_1, fasta_file_2)
    return alignment_score, alignments


class Gotoh:

    def __init__(self, file_substitution_matrix, scores, use_dna):
        self.d_matrix = [[]]
        self.p_matrix = [[]]
        self.q_matrix = [[]]
        self.score_matrix = []
        self.lookup = {}
        self.alpha = scores['alpha']
        self.beta = scores['beta']
        self.use_dna = use_dna
        if not self.use_dna:
            self.set_up_substitution_matrix(file_substitution_matrix)
        # counters for traceback over 3 matrices
        self.i = 0
        self.j = 0
        self.computedAlignment = []
        self.paths = [()]
        self.all_traces = [[]]  # Track the actions diag, left up form pQD to be returned

    def run(self, fasta_file_1, fasta_file_2):
        """Put your code here"""
        # seq1 = self.read_fasta_file(fasta_file_1)
        # seq2 = self.read_fasta_file(fasta_file_2)
        seq1 = 'CG'
        seq2 = 'CCGA'
        self.p_matrix = self.init_matrix_p(seq1, seq2)
        self.q_matrix = self.init_matrix_q(seq1, seq2)
        self.d_matrix = self.init_matrix_d(seq1, seq2, self.alpha, self.beta)
        self.complete_d_p_q_computation(seq1, seq2, self.alpha, self.beta)
        trace_list = self.compute_all_tracebacks(seq1, seq2, self.d_matrix,
                                                 self.p_matrix, self.q_matrix, self.alpha, self.beta)
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

    # get score

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

    def dna_match_mismatch(self, char1, char2):
        """
         Args:
            char1: character from seq1
            char2: character from seq1

        Returns:
            score based on the match/mismatch
        """
        if char1 == char2:
            return SCORES_DNA['match']
        elif char1 != char2:
            return SCORES_DNA['mismatch']

    def get_score(self, d, char1, char2):
        if self.use_dna:
            return d + self.dna_match_mismatch(char1, char2)
        else:
            return d + self.read_substitution_matrix(char1, char2)

    def affine_gap(self, i):
        return self.alpha + i * self.beta  # SCORES_DNA["alpha"] + SCORES_DNA["beta"]*i  # g(k) = -3 - k

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
            matrix_d[i][0] = self.affine_gap(i)
        for j in range(1, m):
            matrix_d[0][j] = self.affine_gap(j)
        return matrix_d

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

    def calculate_p(self, value_d, value_p):
        return max(value_d + self.affine_gap(1), value_p + self.beta)

    def calculate_q(self, value_d, value_q):
        return max(value_d + self.affine_gap(1), value_q + self.beta)

    def calculate_d(self, value_d, value_p, value_q, char1, char2):
        new_d = self.get_score(value_d, char1, char2)
        return max(value_p, max(value_q, new_d))

    def complete_d_p_q_computation(self, seq_1, seq_2, cost_gap_open, cost_gap_extend, substitutions=None):
        """
        Implement the recursive computation of matrices D, P and Q
        """
        n = len(seq_1) + 1
        m = len(seq_2) + 1
        for i in range(1, n):
            for j in range(1, m):
                self.p_matrix[i][j] = self.calculate_p(self.d_matrix[i - 1][j], self.p_matrix[i - 1][j])
                self.q_matrix[i][j] = self.calculate_q(self.d_matrix[i][j - 1], self.q_matrix[i][j - 1])
                self.d_matrix[i][j] = self.calculate_d(self.d_matrix[i - 1][j - 1], self.p_matrix[i][j],
                                                       self.q_matrix[i][j], seq_1[i - 1], seq_2[j - 1])
                self.visualize_matrix(self.p_matrix)
                self.visualize_matrix(self.q_matrix)
                self.visualize_matrix(self.d_matrix)

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
        # using structure [((row,column),"matrix")] for traceback
        self.i = len(d_matrix) - 1
        self.j = len(d_matrix[0]) - 1
        self.trace_ctr = 0  # for copys path
        self.paths = [((self.i, self.j), 'D')]  # track what matrix and row/col for matrix traversal
        self.all_traces = [[]]  # Track the actions diag, left up form pQD to be returned

        while self.i > 0 or self.j > 0:
            if self.paths[self.trace_ctr][1] == 'D':
                # check on main matrix
                self.trace_d(seq1[self.i - 1], seq2[self.j - 1])
            elif self.paths[self.trace_ctr][1] == 'P':
                self.trace_p()
            elif self.paths[self.trace_ctr][1] == 'Q':
                self.trace_q()
            (self.i, self.j) = self.paths[self.trace_ctr][0]

        return self.all_traces  # all_paths

    def trace_d(self, char1, char2):
        # var to track i, j local
        local_i = self.i
        local_j = self.j
        split = False
        if self.j > 0 and self.i > 0:
            if self.d_matrix[self.i][self.j] == \
                    self.d_matrix[self.i - 1][self.j - 1] + self.dna_match_mismatch(char1, char2):
                print(f"came from D matrix diagonal")
                self.all_traces[self.trace_ctr].append('diag_D')
                local_i -= 1
                local_j -= 1
                split = False
            if self.d_matrix[self.i][self.j] == self.q_matrix[self.i][self.j]:  # possible other way
                print(f"SPlit from Q ,go to Q")

            if self.d_matrix[self.i][self.j] == self.p_matrix[self.i][self.j]:  # another way
                print(f"split from P ,go to P")

    def trace_p(self):
        # var to track i, j local
        split = False
        if self.i > 0:
            if self.p_matrix[self.i][self.j] == self.d_matrix[self.i - 1][self.j] + self.affine_gap(1):
                self.all_traces[self.trace_ctr].append('up_D')
                i, j = self.paths[self.trace_ctr][0]
                self.paths[self.trace_ctr] = ((i - 1, j), 'D')  # <--split happened
                split = True
            if self.p_matrix[self.i][self.j] == self.p_matrix[self.i - 1][self.j] + SCORES_DNA['beta']:
                if split:
                    self.all_traces.append(self.all_traces[self.trace_ctr][0:-1])
                    self.all_traces[len(self.all_traces) - 1].append('up_P')
                    self.paths.append(((self.i - 1, self.j), 'P'))
                else:
                    self.all_traces[self.trace_ctr].append('up_P')
                    i, j = self.paths[self.trace_ctr][0]
                    self.paths[self.trace_ctr] = ((i - 1, j), 'P')

    def trace_q(self):
        split = False
        if self.j > 0:
            if self.q_matrix[self.i][self.j] == self.d_matrix[self.i][self.j - 1] + self.affine_gap(1):
                self.all_traces[self.trace_ctr].append('left_D')
                i, j = self.paths[self.tracebackStackIndex][0]
                self.paths[self.trace_ctr] = ((i, j - 1), 'D')  # <--split happened
                split = True

            if self.q_matrix[self.i][self.j] == self.q_matrix[self.i][self.j - 1] + SCORES_DNA['beta']:
                if split:
                    self.all_traces.append(self.all_traces[self.trace_ctr][0:-1])
                    self.all_traces[len(self.all_traces) - 1].append('left_Q')
                    self.paths.append(((self.i, self.j - 1), 'Q'))
                else:
                    self.all_traces[self.trace_ctr].append('left_Q')
                    self.paths[self.trace_ctr] = ((self.i, self.j - 1), 'Q')

    def alignment(self, traceback, seq1, seq2):
        """
        Implement creation of the alignment with given traceback path and sequences1 and 2
        """
        i = len(seq1) - 1
        j = len(seq2) - 1
        seq1_align = ""
        seq2_align = ""
        for step in traceback:
            if step == "diag_D":
                seq1_align += seq1[i]
                seq2_align += seq2[j]
                i -= 1
                j -= 1

            elif step == "up_D" or step == 'up_P':
                seq1_align += seq1[i]
                seq2_align += '-'
                i -= 1

            elif step == "left_D" or step == 'left_Q':
                seq1_align += '-'
                seq2_align += seq2[j]
                j -= 1

        while j >= 0:
            seq1_align += '-'
            seq2_align += seq2[j]
            j -= 1

        while i >= 0:
            seq2_align += '-'
            seq1_align += seq1[i]
            i -= 1
        return seq1_align[::-1], seq2_align[::-1]

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
    # DNA
    gotoh(fasta1, fasta2, SCORES_DNA, True)
    # protien pam
    # gotoh(fasta1, fasta2, SCORES_PRO, False, pam_file)
    # protien blosum
    # gotoh(fasta1, fasta2, SCORES_PRO, False, blosum_file)
    # Tool to test o/p http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh
