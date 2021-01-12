
def gotoh(fasta_file_1, fasta_file_2, cost_gap_open, file_substitution_matrix=None):
    """Put your code here"""
    algo = Gotoh()
    alignment_score, alignments = algo.run(fasta_file_1, fasta_file_2, cost_gap_open, file_substitution_matrix)
    return alignment_score, alignments


class Gotoh:

    def run(self, fasta_file_1, fasta_file_2, cost_gap_open, file_substitution_matrix=None):
        """Put your code here"""
        seq1 = self.read_fasta_file(fasta_file_1)
        seq2 = self.read_fasta_file(fasta_file_2)
        sub_matrix = self.read_substitution_matrix(file_substitution_matrix)
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

    def read_substitution_matrix(self, file_substitution_matrix):
        """
        Implement reading the scores file.
        It can be stored as a dictionary of example:
        scores[("A", "R")] = -1

        Args:
            file_substitution_matrix: location of substituition matrix

        Returns:
            substituition matrix

        """
        for line in open(file_substitution_matrix):
            li = line.strip()
            if not li.startswith("#"):
                print(line.rstrip())

        return 0  # scores

    def init_matrix_d(self, seq_1, seq_2, cost_gap_open, cost_gap_extend):
        """
        Implement initialization of the matrix D
        """
        return []  # matrix_d

    def init_matrix_p(self, seq_1, seq_2):
        """
        Implement initialization of the matrix P
        """
        return []  # matrix_p

    def init_matrix_q(self, seq_1, seq_2):
        """
        Implement initialization of the matrix Q
        """
        return []  # matrix_q

    def complete_d_p_q_computation(self, seq_1, seq_2, cost_gap_open, cost_gap_extend, substitutions=None):
        """
        Implement the recursive computation of matrices D, P and Q
        """

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
        return []  # all_paths

    def find_all_previous(self, cell, seq1, seq2, d_matrix, p_matrix, q_matrix,
                          cost_gap_open, cost_gap_extend, substitution=None):
        parent_cells = []
        """
        Implement a search for all possible previous cells.
        """

        return parent_cells

    def check_complete(self, path):
        """
        Implement a function which checks if the traceback path is complete.
        """

    def alignment(self, traceback_path, seq1, seq2):
        """
        Implement creation of the alignment with given traceback path and sequences1 and 2
        """
        return '-', '-'  # alignment_seq1, alignment_seq2

    # HELPER FUNCTIONS FOR SELF CHECK

    def visualize_matrix(self, matrix):
        """
        Implement the visualization of a matrix.
        Can be used for self check
        """

    def score_of_alignment(self, align_seq1, align_seq2, cost_gap_open,
                           cost_gap_extension, substitutions=None):
        """
        A nice helper function which computes the score of the given alignment.
        This is only used for the self check.
        Input example:
        --CG
        AACA
        """
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


