import math


def check_input(alignments):
    """Checks if the input is correct. Raises a ValueError if it is not."""
    # the input has to be a dictionary
    if not isinstance(alignments, dict):
        raise ValueError("malformed input")

    for key in alignments:
        # key must be a tuple with 2 numbers
        if not isinstance(key, tuple) or len(key) != 2:
            raise ValueError("malformed input")
        i = key[0]
        j = key[1]

        # the numbers have to be ints (True and False dont count!)
        if not is_good_number(i) or not is_good_number(j):
            raise ValueError("malformed input")
        if i >= j:
            raise ValueError("malformed input")

        # the value should be 2 strings
        value = alignments[key]
        if not isinstance(value, tuple) or len(value) != 2:
            raise ValueError("malformed input")
        seq1 = value[0]
        seq2 = value[1]
        if not isinstance(seq1, str) or not isinstance(seq2, str):
            raise ValueError("malformed input")
        if len(seq1) != len(seq2):
            raise ValueError("malformed input")

        # only A, C, G, T and - are allowed
        bad_letters = [c for c in (seq1 + seq2).upper() if c not in "ACGT-"]
        if len(bad_letters) > 0:
            raise ValueError("malformed input")


def is_good_number(x):
    """Returns True if x is an int that is 0 or bigger (but not True/False)."""
    if isinstance(x, bool):
        return False
    if not isinstance(x, int):
        return False
    return x >= 0


def get_number_of_sequences(alignments):
    """Returns how many sequences there are. Raises a ValueError if
    some pairs are missing."""
    if len(alignments) == 0:
        return 1  # only one sequence, so no pairs

    n = max([key[1] for key in alignments]) + 1

    # there should be n*(n-1)/2 pairs, otherwise something is missing
    if len(alignments) != n * (n - 1) // 2:
        raise ValueError("malformed input")
    return n


def calculate_distance(seq1, seq2):
    """Calculates the Jukes-Cantor distance between two aligned sequences.
    Positions with a gap are skipped. Returns 30.0 if the distance
    can not be calculated."""
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    # only use the positions where there is no gap
    positions = [k for k in range(len(seq1)) if seq1[k] != "-" and seq2[k] != "-"]
    if len(positions) == 0:
        return 30.0  # nothing to compare

    differences = [k for k in positions if seq1[k] != seq2[k]]
    p = len(differences) / len(positions)

    if p >= 0.75:
        return 30.0  # formula doesnt work here (log of negative number)
    return -0.75 * math.log(1 - (4/3) * p)


def ComputeDistMatrix(alignments):
    """Returns a matrix with the evolutionary distances between all
    sequences. Raises a ValueError if the input is wrong."""
    check_input(alignments)
    n = get_number_of_sequences(alignments)

    # make a matrix full of zeros
    matrix = [[0.0 for b in range(n)] for a in range(n)]

    # fill in the distance for every pair
    for key in alignments:
        i = key[0]
        j = key[1]
        d = calculate_distance(alignments[key][0], alignments[key][1])
        matrix[i][j] = d
        matrix[j][i] = d

    return matrix