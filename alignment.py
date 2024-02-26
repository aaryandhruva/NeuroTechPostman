import sys

#Match/Mismatch: The score is incremented if the characters at the corresponding positions in both sequences match.
#Gap (Insertion or Deletion): The score is incremented when introducing a gap in one of the sequences to align them.
#Extension of an existing gap: Sometimes, it's allowed to extend a gap without incrementing the score.
def read_sequence(file_path):
    with open(file_path, 'r') as file:
        return file.readline().strip()

def pairwise_alignment(seq1, seq2):
    # initialising the alignment score matrix
    score_matrix = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    # fillung the score matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            match = score_matrix[i - 1][j - 1] + (1 if seq1[i - 1] == seq2[j - 1] else 0)
            delete = score_matrix[i - 1][j]
            insert = score_matrix[i][j - 1]
            score_matrix[i][j] = max(match, delete, insert)

    # traceback to find the alignment
    alignment_seq1 = ''
    alignment_seq2 = ''
    i, j = len(seq1), len(seq2)
    while i > 0 and j > 0:
        if seq1[i - 1] == seq2[j - 1]:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            i -= 1
            j -= 1
        elif score_matrix[i][j] == score_matrix[i - 1][j] + 1:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = '-' + alignment_seq2
            i -= 1
        else:
            alignment_seq1 = '-' + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            j -= 1

    #  if one sequence is longer than the other
    while i > 0:
        alignment_seq1 = seq1[i - 1] + alignment_seq1
        alignment_seq2 = '-' + alignment_seq2
        i -= 1
    while j > 0:
        alignment_seq1 = '-' + alignment_seq1
        alignment_seq2 = seq2[j - 1] + alignment_seq2
        j -= 1

    return alignment_seq1, alignment_seq2, score_matrix[-1][-1]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("pls enter: python alignment.py <file1> <file2>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    seq1 = read_sequence(file1_path)
    seq2 = read_sequence(file2_path)

    alignment_seq1, alignment_seq2, score = pairwise_alignment(seq1, seq2)

    print("Sequence 1:", alignment_seq1)
    print("Sequence 2:", alignment_seq2)
    print("Alignment Score:", score)
