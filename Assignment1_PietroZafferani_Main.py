''' Designing an algorithm for The Optimal Sequences Alignment '''
'''This module implements the functions for the backtracking and for printing the alignment'''

from Assignment1_PietroZafferani_Matrixes import *

'''In the final alignment the mismatches will be displayed in lowercase and the matches in uppercase '''

def mismatch(a: str, b: str) -> tuple:
    if a != b:
        return (a.lower(), b.lower())   # mismatching characters
    else:
        return (a.upper(), b.upper())   # matching characters


'''Adds insertions in one of the 2 sequences aligned'''

def Indels() -> str:
    return '-'      # insertion symbol


'''This function backtracks in the Matrix containing the optimal path for the GLOBAL alignment, every case is 
    taken into account and the relative characters are added to the newly created aligned sequences. this
    implementation yields the reverse of the actual sequences and they can be easily printed correctly by another function.
    The recursion takes linear time to reach the start of the sequences.'''
# M[0][coordinate2]--> are the labels of first sequence (x-axis)
# M[coordinate1][0] --> are the labels for second sequence(y-axis)

def backtracking(M: list, coordinate1: int, coordinate2: int, x='', y='') -> tuple:  # c1 = rows # c2 = columns

    if coordinate1 == 1 and coordinate2 == 1:   # base case of recursion
        return (x, y)  # first position is seq1(x-axis) #second position is seq2(y-axis)

    elif coordinate2 == 1:  # we reached the extra column so we can only go upwards
        x += Indels()
        y += M[coordinate1][0].upper()
        return backtracking(M, coordinate1 - 1, coordinate2, x, y)

    elif coordinate1 == 1:      # we reached the extra row so we can only go leftwards
        x += M[0][coordinate2].upper()
        y += Indels()
        return backtracking(M, coordinate1, coordinate2 - 1, x, y)

    elif M[coordinate1][coordinate2] == 'DIAG':     # calling the recursion in diagonal
        evaluate = mismatch(M[0][coordinate2], M[coordinate1][0])
        x += str(evaluate[0])
        y += str(evaluate[1])
        return backtracking(M, coordinate1 - 1, coordinate2 - 1, x, y)

    elif M[coordinate1][coordinate2] == 'UP':       # calling the recursion in the cell above
        x += Indels()
        y += M[coordinate1][0].upper()
        return backtracking(M, coordinate1 - 1, coordinate2, x, y)

    elif M[coordinate1][coordinate2] == 'LEFT':     # calling the recursion in the previous cell on the left
        x += M[0][coordinate2].upper()
        y += Indels()
        return backtracking(M, coordinate1, coordinate2 - 1, x, y)


'''Defines the parameters to be used in the GLOBAL optimal alignment and invokes the traceback'''

def Global(seq1: str, seq2: str, gap: int, m: int) -> tuple:
    Matrixes = fillMatrix(seq1, seq2, gap, gap, m)  # in the global alignment the gap coincides with the extra gap
    Backtrack = Matrixes[2]     # matrix storing all the moves

    return backtracking(Backtrack, len(Backtrack) - 1, len(Backtrack[0]) - 1)
    # the two coordinates correspond always to the last cell in the Matrix (down-right corner)


'''Same backtracking function but implemented for the LOCAL optimal alignment'''
# M[0][coordinate2]--> are the labels of first sequence (x-axis)
# M[coordinate1][0] --> are the labels for second sequence(y-axis)

def Local_backtracking(Mparal: list, Mnumbers: list, coordinate1: int, coordinate2: int, x='', y='') -> tuple:  # coordinate1=rows #coordinate2=columns

    if coordinate1 == 1 and coordinate2 == 1:  # we reached the first cell of the Matrix so we simply return the result
        return ('<' + x + '>', '<' + y + '>')   # base case of recursion

    elif Mnumbers[coordinate1][coordinate2] == 0:  # other base case of recursion
        if coordinate1 == 1:            # sub-case when the extra row is reached
            x += Mparal[0][coordinate2].upper()
            y += Indels()
        elif coordinate2 == 1:          # sub-case when the extra column is reached
            x += Indels()
            y += Mparal[coordinate1][0].upper()
        else:                           # sub-case when the null-cell is the 'middle' of the matrix
            evaluate = mismatch(Mparal[0][coordinate2], Mparal[coordinate1][0])
            x += str(evaluate[0])
            y += str(evaluate[1])

            # consider that everything must be reversed
        return ('<' + x + '>', '<' + y + '>')  # we reached the end of the recursion

    elif Mparal[coordinate1][coordinate2] == 'DIAG':  # recursion in diagonal
        evaluate = mismatch(Mparal[0][coordinate2], Mparal[coordinate1][0])
        x += str(evaluate[0])
        y += str(evaluate[1])
        return Local_backtracking(Mparal, Mnumbers, coordinate1 - 1, coordinate2 - 1, x, y)

    elif Mparal[coordinate1][coordinate2] == 'UP':  # recursion in vertical
        x += Indels()
        y += Mparal[coordinate1][0].upper()
        return Local_backtracking(Mparal, Mnumbers, coordinate1 - 1, coordinate2, x, y)

    elif Mparal[coordinate1][coordinate2] == 'LEFT':  # recursion in horizontal
        x += M[0][coordinate2].upper()
        y += Indels()
        return Local_backtracking(Mparal, Mnumbers, coordinate1, coordinate2 - 1, x, y)


'''Defines the parameters to be used in the LOCAL optimal alignment and invokes the traceback'''

def Local(seq1: str, seq2: str, gap: int, m: int) -> tuple:
    Matrixes = fillMatrix(seq1, seq2, gap, 0, m)
    Mnumbers = Matrixes[0]      # matrix containing all the values of the alignments
    coordinates = Matrixes[1]   # tuple with the coordinates of the highest-value cell
    Mparal = Matrixes[2]     # matrix storing the moves performed to fill the cells

    return Local_backtracking(Mparal, Mnumbers, coordinates[0], coordinates[1])


'''This function iterates over the tuple containing the two sequences(reversed) and prints them one
    on top of each other'''

def PrintAlignment(aligned_tuple: tuple) -> print:
    for sequence in aligned_tuple:
        print(sequence[::-1])   # the strings were created in the opposite sense


'''This is the only function that the user must be call to interact with the system'''

def MAIN()-> print:
    seq1 = input('Please type in the first sequence you would like to align: '+'\n')
    seq2 = input('Please type in the second sequence you would like to align: '+'\n')
    type = input('What kind of alignment would you like to perform:(select GLOBAL/LOCAL) '+'\n')
    gap = input('Please select the score assigned to gaps: (suggested -2)'+'\n')
    score= input('Please select the score assigned to a match of characters: (suggested 1)'+'\n')

    print('\n' + type.upper() + ' alignment:' + '\n')

    if type.upper() == 'GLOBAL':
        result = Global( str(seq1), str(seq2), int(gap), int(score))
        PrintAlignment(result)

    else:
        result = Local( str(seq1), str(seq2), int(gap), int(score))
        PrintAlignment(result)


''''This section of the file is designed to be the only interacting part with the user'''

if __name__ == '__main__':

    '''Below are given some biological sequences that can be use to test the algorithms
     by simply coping and pasting them properly.'''

    #> CAA25109.1 myoglobin[Homo sapiens]
    'MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG'

    #> NP_999401.1 myoglobin[Sus scrofa]
    'MGLSDGEWQLVLNVWGKVEADVAGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGNTVLTALGGILKKKGHHEAELTPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMSKALELFRNDMAAKYKELGFQG'


    #NB: the quotation marks '' are not needed when typing the input;
       # only letters are valid as the argument of the aligment

    MAIN()