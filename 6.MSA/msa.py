import re
import numpy as np
import operator
from new_gotoh import Gotoh

dnaA = 'AAAAATTTACATTTAA'

dnaB = 'ATTTTACATTTGGG'

dnaC = 'ATTTTACATTTGGGG'

sample1 = "CG"
sample2 = "CCGA"
sample3 = 'CCAG'


def create_group(old_align, new_align):
    new_grp =[]
    for al in old_align:
        new_grp.append(re.sub('-', '#', al))
    for al in new_align:
        new_grp.append(re.sub('-', '#', al))
    return new_grp


if __name__ == '__main__':
    algo = Gotoh()
    # pick the shortest as the initial first group.
    alignments = [dnaA,dnaB,dnaC]
    shortest = min(alignments, key=len)
    alignments = np.setdiff1d(alignments, shortest)

    score_3 = []
    for al in alignments:
        scores, alignment = algo.run(shortest, al)
        align = {
            'score': scores,
            'alignments': alignment,
            'seq1': shortest,
            'seq2': al
        }
        score_3.append(align)
        algo.reset()

     # given the alignments get the one with the max scored alignments with all info
    max_alignmentgrp = sorted(score_3, key=lambda x: (x['score'], x['alignments'], x['seq1'], x['seq2']) , reverse = True)[0]

    add_seq1 = max_alignmentgrp['alignments'][0][0]
    add_seq2 = max_alignmentgrp['alignments'][0][1]
    group1 = create_group([add_seq1],[add_seq2])
    score_3 = []
    # compare b/w A and members of the group
    for mem in group1:
        scores, alignment = algo.run(mem, dnaA)
        align = {
            'score': scores,
            'alignments': alignment,
            'seq1': dnaA,
            'seq2': mem
        }
        score_3.append(align)
        algo.reset()

    max_alignmentgrp = sorted(score_3, key=lambda x: (x['score'], x['alignments'], x['seq1'], x['seq2']) ,reverse = True)[0]
    new_al = np.setdiff1d(group1, [max_alignmentgrp['seq1'], max_alignmentgrp['seq2']])
    # convert the tuple to list
    new_group = create_group(list(max_alignmentgrp['alignments'][0]), new_al)

    print(f"final grouping{new_group}")

