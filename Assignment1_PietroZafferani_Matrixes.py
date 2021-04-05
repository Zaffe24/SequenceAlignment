''' Implementing The Needleman & Wunsh and The Smith & Waterman algorithms for The Optimal Sequences Alignment '''
    # by Pietro Zafferani, Genomics

'''This module implements functions to create and to handle alignment matrix and the backtracking one'''

'''This function creates the skeleton of the both matrixes, only the columns and rows indexes are present.'''

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
    the extra gaps are all zeros for the local one (extra_gap=0);
    the extra gaps sum on top of each other for the global one (extra_gap=integer)'''

def typeMatrix( seq1: str, seq2: str, extra_gap: int) -> list:
    M = emptyMatrix(seq1, seq2)  # create the skeleton
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

def match(i, j, n: int) -> int:
    return n if i == j else -n


'''This function is used to assign the best moove among the three possibilities'''

def get_max(diag: int, up: int, left: int) -> str:
    if diag >= up and diag >= left:  # the diagonal moove is favoured in case of equality since a mismatch is lighter than a gap
        return 'DIAG'
    elif up >= diag and up >= left:
        return 'UP'
    elif left >= diag and left >= up:
        return 'LEFT'


'''This function returns a tuple containing the matrix filled with all the scores,
    the coordinates of the cell with highest value (needed for the local implementation)
    and the matrix storing all the mooves to reach a certain cell'''

def fillMatrix(seq1: str, seq2: str, gap: int, extra_gap: int, m: int) -> tuple:
    coordinates = (1, 1)  # default value in case there is no a positive cell's value in the whole matrix
    Matrix = typeMatrix(seq1, seq2, extra_gap)       # numerical matrix
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

            cell_value = max(diag, up, left)    # in case of a tie for the highest value the first one is chosen
            paralcel = get_max(diag, up, left)  # this cell keeps track of the moves in the backtracking matrix

            Matrix[i].append(cell_value)     # append the new cell to the current row
            ParalMatrix[i].append(paralcel)  # storing the move to reach the current cell

            if cell_value >= best_value:    # feature needed to implement local optimal alignment, it saves the coordinates of the highest-value cell
                best_value = cell_value
                coordinates = (i, j)

            j += 1  # move to the next cell on the same row
        i += 1      # move to the next row

    return (Matrix, coordinates, ParalMatrix)  # this tuple contains all the info for the backtracking function


'''The following functions allow to visualize graphically the matrixes,
    not actually functional to the algorithm but it's helpful in the testing'''

'''Prints the matrix containing all the values'''

def SHOW(M) -> print:  # .fillMatrix() as parameter
    for i in M[0]:
        print(i)
    print('highest cell coordinates :' + str(M[1]))


'''Prints the matrix storing all the moves'''

def SHOWparallel(M) -> print:  # .fillMatrix() as parameter
    for i in M[2]:
        print(i)


if __name__ == '__main__':

    #seq1 = 'xxxxx'
    #seq2 = 'vvvvvvvvvvvvvvvv'
    #gap = -2 #negative value
    #extra_gap = 0 #negative value or 0 according to the type of alignment
    #score = 1 #absolute value of match and mismatch

    '''call these functions to print the matrixes'''

    #SHOWparallel(fillMatrix(seq1,seq2,gap,extra_gap,score))
    #print()
    #SHOW(fillMatrix(seq1,seq2,gap,extra_gap,score))
