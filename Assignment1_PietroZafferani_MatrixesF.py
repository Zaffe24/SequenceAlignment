'''Designing an algorithm for The Optimal Sequences Alignment'''
'''This module implements functions to create and to handle alignment matrix and the backtracking one'''

'''This function creates the skeleton of the alignment matrix, only the columns and rows indexes are present.'''

def emptyMatrix(seq1: str, seq2: str) -> list:
    Matrix = [['/', '/']]  # extra null column
    for column in seq1:
        Matrix[0].append(column)
    Matrix.append(['/'])  # extra null row
    for row0 in seq2:
        Row = [row0]
        Matrix.append(Row)

    return Matrix


'''This function defines the type of alignment:
    the extra gaps are all zeros for the local one (extra_gap=0), set by default;
    the extra gaps sum on top of each other for the global one (extra_gap=value)'''

def typeMatrix(seq1: str, seq2: str, extra_gap: int) -> list:
    M = emptyMatrix(seq1, seq2)
    start = 0
    for column in range(len(seq1) + 1):
        M[1].append(start)
        start += extra_gap  # filling the extra null row

    start = extra_gap
    for row in M[2:]:
        row.append(start)
        start += extra_gap  # filling the extra null column

    return M  # the type of matrix depends on the value of the extra_gap:
    # 0 --> for local alignment ; (-)integer --> for global alingment


'''function used to determine the score of matches/mismatches'''

def match(i, j, n: int):
    return n if i == j else -n


'''This function is used to assign the best moove among the three possibilities'''

def get_max(diag: int, up: int, left: int) -> str:
    if diag >= up and diag >= left:  # the diagonal moove is favoured in case of equality since a mismatch is lighter than a gap
        return 'DIAG'
    elif up >= diag and up >= left:
        return 'UP'
    elif left >= diag and left >= up:
        return 'LEFT'


def fillMatrix(seq1: str, seq2: str, gap: int, extra_gap: int, m: int) -> tuple:
    coordinates = (0, 0)  # default value
    Matrix = typeMatrix(seq1, seq2, extra_gap)      # numerical matrix
    ParalMatrix = typeMatrix(seq1, seq2, extra_gap)  # backtracking matrix
    i = 2
    best_value = 0
    for row in Matrix[i:]:  # start iterating from the first non-extra row until we reach the last character of seq2
        j = 2
        for cell in Matrix[0][j:]:   # adding new cells corresponding to the number of characters of seq1

            # thre possible mooves to reach the cell
            diag = Matrix[i - 1][j - 1] + match(row[0], Matrix[0][j], m)
            up = Matrix[i - 1][j] + gap
            left = row[j - 1] + gap

            cell_value = max(diag, up, left)  # in case of a tie for the highest value the first one is chosen
            paralcel = get_max(diag, up, left)  # this cell keeps track of the moves in the backtracking matrix

            Matrix[i].append(cell_value)    # append the new cell to the current row
            ParalMatrix[i].append(paralcel) # storing the moove to reach the current cell

            if cell_value >= best_value:  # feature needed to implement local optimal alignment, we save the coordinates of the highest cell value
                best_value = cell_value
                coordinates = (i, j)

            j += 1  # moove to the next cell on the same row
        i += 1      # moove to the next row

    return (Matrix, coordinates, ParalMatrix) # this tuple contains all the info for the backtacking function


'''The following functions allow to visualize graphically the matrixes,
    not actually functional to the algorithm but it's helpful in the testing'''

def SHOW(M):
    for i in M[0]:
        print(i)
    print('highest cell coordinates :' + str(M[1]))


def SHOWparal(M):
    for i in M[2]:
        print(i)


'''This section of the file interacts with the user'''
if __name__ == '__main__':
    '''Below are reported some biological sequences suitable for testing the algorithms'''

    # >XM_019013136.2 PREDICTED: Gorilla gorilla gorilla tumor protein p53 (TP53), transcript variant X1, mRNA (1:80 nts)
    GorillaP53 = 'GAATTAAAATAGGATGACTTAAAGTCTGCACGGGAAGGAGCCTACCCCCATGTTCCTGGCTAGCCAAGGAACCACCAGTT'

    # >NM_000546.6 Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA(1:80 nts)
    HumnaP53 = 'CTCAAAAGTCTAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTG'

    A = emptyMatrix('ABCD', 'ABCD')
    a = typeMatrix('ABCDE', 'ABCDE', -2)
    fullM = fillMatrix('t', 'AA', -2, -2, 1)
    SHOW(fullM)
    # SHOWparal(fullM)
