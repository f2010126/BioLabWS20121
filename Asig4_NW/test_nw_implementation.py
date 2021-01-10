import nw_implementation as nw
import numpy as np

seq1 = 'AAA'
seq2 = 'GGG'
def test_calculate_matrix():
    assert True


def test_init_matrix():
    test = [[ 0., -1., -2., -3.],
                   [-1.,  0.,  0.,  0.],
                   [-2.,  0.,  0.,  0.],
                   [-3.,  0.,  0.,  0.]]
    algo = nw.NeedlemanWunsch()
    init = algo.init_matrix(seq1,seq2, -1)
    assert np.array_equal(init, test), "Initialisation error"


def test_read_fasta_file():
    assert True
