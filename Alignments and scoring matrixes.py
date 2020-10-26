"""
Code for analysing pairwise sequence alignments for DNA, RNA and protein folds with the use of global and local alignments
"""
import math
import random
import matplotlib.pyplot as plt


# functions

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    outputs a dictionary of dictioraries representing a scoring matrix
    """
    score = {}
    characters = set(alphabet)
    characters.add('-')
    for letter in characters:
        score[letter] = {}
        for alignment in characters:
            if letter == '-' or alignment == '-':
                score[letter][alignment] = dash_score
            elif letter == alignment:
                score[letter][alignment] = diag_score
            else:
                score[letter][alignment] = off_diag_score
    return score

def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    returns the alignemnt matrix for seq_x and seq_y, if global_flag is true we perform a global alignment or, otherwise, a local alignemnt matrix
    """
    len_x = len(seq_x) + 1
    len_y = len(seq_y) + 1
    dynamic_table = [[0 for col in range(len_y)] for row in range(len_x)]
    dynamic_table[0][0] = 0

    for index in range(1, len_x):
        dynamic_table[index][0] = dynamic_table[index - 1][0] + \
            scoring_matrix[seq_x[index - 1]]['-']
        if not global_flag and dynamic_table[index][0] < 0:
            dynamic_table[index][0] = 0

    for index in range(1, len_y):
        dynamic_table[0][index] = dynamic_table[0][index - 1] + \
            scoring_matrix['-'][seq_y[index - 1]]
        if not global_flag and dynamic_table[0][index] < 0:
            dynamic_table[0][index] = 0

    for row in range(1, len_x):
        for col in range(1, len_y):
            alignments = (dynamic_table[row - 1][col - 1] +
                          scoring_matrix[seq_x[row - 1]][seq_y[col - 1]],
                          dynamic_table[row - 1][col] +
                          scoring_matrix[seq_x[row - 1]]['-'],
                          dynamic_table[row][col - 1] +
                          scoring_matrix['-'][seq_y[col - 1]])
            dynamic_table[row][col] = max(alignments)
            if not global_flag and dynamic_table[row][col] < 0:
                dynamic_table[row][col] = 0
    return dynamic_table

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    returns a tuple with the optimal score of the alignment and the alignments for sequence x and sequence y
    """
    if seq_x == '' or seq_y == '':
        return (0, seq_x, seq_y)

    len_x = len(seq_x)
    len_y = len(seq_y)
    alignment_x = ''
    alignment_y = ''
    score = 0

    while len_x != 0 and len_y != 0:
        if alignment_matrix[len_x][len_y] == alignment_matrix[len_x - 1][len_y - 1] + scoring_matrix[seq_x[len_x - 1]][seq_y[len_y - 1]]:
            alignment_x = seq_x[len_x - 1] + alignment_x
            alignment_y = seq_y[len_y - 1] + alignment_y
            score += scoring_matrix[seq_x[len_x - 1]][seq_y[len_y - 1]]
            len_x -= 1
            len_y -= 1
        else:
            if alignment_matrix[len_x][len_y] == alignment_matrix[len_x - 1][len_y] + scoring_matrix[seq_x[len_x - 1]]['-']:
                alignment_x = seq_x[len_x - 1] + alignment_x
                alignment_y = '-' + alignment_y
                score += scoring_matrix[seq_x[len_x - 1]]['-']
                len_x -= 1
            else:
                alignment_x = '-' + alignment_x
                alignment_y = seq_y[len_y - 1] + alignment_y
                score += scoring_matrix['-'][seq_y[len_y - 1]]
                len_y -= 1

    while len_x != 0:
        alignment_x = seq_x[len_x - 1] + alignment_x
        alignment_y = '-' + alignment_y
        score += scoring_matrix[seq_x[len_x - 1]]['-']
        len_x -= 1

    while len_y != 0:
        alignment_x = '-' + alignment_x
        alignment_y = seq_y[len_y - 1] + alignment_y
        score += scoring_matrix['-'][seq_y[len_y - 1]]
        len_y -= 1
    return (score, alignment_x, alignment_y)


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    returns a tuple with the optimal score of the alignment and the alignments for sequence x and sequence y
    """
    if seq_x == '' or seq_y == '':
        return (0, seq_x, seq_y)

    score, (idx1, idx2) = max((score, (idx1, idx2))
                              for idx1, row in enumerate(alignment_matrix) for idx2, score in enumerate(row))
    alignment_x = ''
    alignment_y = ''
    score = 0

    while idx1 > 0 and idx2 > 0:
        if alignment_matrix[idx1][idx2] == alignment_matrix[idx1 - 1][idx2 - 1] + scoring_matrix[seq_x[idx1 - 1]][seq_y[idx2 - 1]]:
            alignment_x = seq_x[idx1 - 1] + alignment_x
            alignment_y = seq_y[idx2 - 1] + alignment_y
            score += scoring_matrix[seq_x[idx1 - 1]][seq_y[idx2 - 1]]
            idx1 -= 1
            idx2 -= 1
            if alignment_matrix[idx1][idx2] == 0:
                break
        else:
            if alignment_matrix[idx1][idx2] == alignment_matrix[idx1 - 1][idx2] + scoring_matrix[seq_x[idx1 - 1]]['-']:
                alignment_x = seq_x[idx1 - 1] + alignment_x
                alignment_y = '-' + alignment_y
                score += scoring_matrix[seq_x[idx1 - 1]]['-']
                idx1 -= 1
                if alignment_matrix[idx1][idx2] == 0:
                    break
            else:
                alignment_x = '-' + alignment_x
                alignment_y = seq_y[idx2 - 1] + alignment_y
                score += scoring_matrix['-'][seq_y[idx2 - 1]]
                idx2 -= 1
                if alignment_matrix[idx1][idx2] == 0:
                    break

    return (score, alignment_x, alignment_y)

# data functions


def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = open(filename, 'r')
    ykeys = scoring_file.readline()
    ykeychars = ykeys.split()
    for line in scoring_file.readlines():
        vals = line.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict


def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = open(filename)
    protein_seq = protein_file.read()
    protein_seq = protein_seq.rstrip()
    return protein_seq


def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = open(filename)

    # read in files as string
    words = word_file.read()

    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print("\nLoaded a dictionary with", len(word_list), "words.")
    return word_list


def load_score_matrix(filename):
    file = open(filename)
    return file

# helper functions for analysis questions


def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    """
    Function for returning a dictionary that represents un-normalized
    distribution of num_trial possible sequences of aminoacids
    """
    original_seq_y = list(seq_y)
    scoring_distribution = {}
    for trial in range(1, num_trials + 1):
        random.shuffle(original_seq_y)
        rand_y = ''.join(original_seq_y)
        alignment_matrix = compute_alignment_matrix(seq_x, rand_y,
                                                    scoring_matrix, False)
        local_alignment = compute_local_alignment(seq_x, rand_y,
                                                  scoring_matrix,
                                                  alignment_matrix)
        if local_alignment[0] not in scoring_distribution.keys():
            scoring_distribution[local_alignment[0]] = 0
        scoring_distribution[local_alignment[0]] += 1
    return scoring_distribution


def check_spelling(checked_word, dist, word_list):
    """
    Iterates through word list and returns a set of all words within edit 
    distance (dist) of the string checked word
    """
    similar_words = set([])
    len1 = len(checked_word)
    letters = [chr(letter) for letter in range(97, 123)]
    alphabet = ''.join(letters)
    scoring_matrix = build_scoring_matrix(alphabet, 2, 1, 0)
    for word in word_list:
        alignment_matrix = compute_alignment_matrix(checked_word, word,
                                                    scoring_matrix, True)
        alignment = compute_global_alignment(checked_word, word,
                                             scoring_matrix,
                                             alignment_matrix)
        if len1 + len(word) - alignment[0] <= dist:
            similar_words.add(word)
    return similar_words

# Analysis questions answered


scoring_matrix = read_scoring_matrix("Data//PAM50_scoring_matrix.txt")
human_protein = read_protein("Data//Human_Eyeless_Protein.txt")
fruitfly_protein = read_protein("Data//Fruitfly_Eyeless_Protein.txt")
pax_domain_consensus = read_protein("Data//Consensus_PAX_Domain.txt")

hf_local_alignment_matrix = compute_alignment_matrix(human_protein,
                                                     fruitfly_protein,
                                                     scoring_matrix, False)

hf_local_alignment = compute_local_alignment(human_protein, fruitfly_protein,
                                             scoring_matrix,
                                             hf_local_alignment_matrix)

print("Question 1:\nScore =", hf_local_alignment[0], "and human/fruitfly",
      "sequences are, respectively:\n", hf_local_alignment[1], "\n",
      hf_local_alignment[2])

human_alignment = hf_local_alignment[1]
fruitfly_alignment = hf_local_alignment[2]
human_alignment = human_alignment.replace("-", "")
fruitfly_alignment = fruitfly_alignment.replace("-", "")

hc_global_alignment_matrix = compute_alignment_matrix(human_alignment,
                                                      pax_domain_consensus,
                                                      scoring_matrix, True)

fc_global_alignment_matrix = compute_alignment_matrix(fruitfly_alignment,
                                                      pax_domain_consensus,
                                                      scoring_matrix, True)

hc_global_alignment = compute_global_alignment(human_alignment,
                                               pax_domain_consensus,
                                               scoring_matrix,
                                               hc_global_alignment_matrix)

fc_global_alignment = compute_global_alignment(fruitfly_alignment,
                                               pax_domain_consensus,
                                               scoring_matrix,
                                               fc_global_alignment_matrix)

hc_consensus = [0, 0, 0]
fc_consensus = [0, 0, 0]

for letter1, letter2 in zip(hc_global_alignment[1], hc_global_alignment[2]):
    hc_consensus[0] += 1
    if letter1 == letter2:
        hc_consensus[1] += 1
hc_consensus[2] = hc_consensus[1]/hc_consensus[0]

for letter1, letter2 in zip(fc_global_alignment[1], fc_global_alignment[2]):
    fc_consensus[0] += 1
    if letter1 == letter2:
        fc_consensus[1] += 1
fc_consensus[2] = fc_consensus[1]/fc_consensus[0]

print("\nQuestion 2:\nHuman similarity with PAX domain:",
      str(hc_consensus[2] * 100)
      + "%\nFruitfly similarity with PAX domain:",
      str(fc_consensus[2] * 100) + "%")

print("\nQuestion 3: The level of similarity exhibited by the answers 1 and 2",
      "could have", "been due", "to chance, however that is unlikely to happen.",
      "Due to the vast spectrum of possible combinations of aminoacid bases,",
      "a match that close is statistically improbable.\n")

iterations = 100  # change number of iterations
scoring_distribution = generate_null_distribution(human_protein,
                                                  fruitfly_protein,
                                                  scoring_matrix, iterations)

print("Question 4: score distribution finished calculating", end="")
plt.figure(dpi=1500)
plt.bar(scoring_distribution.keys(), scoring_distribution.values())
plt.title("Distribution of scores of random protein sequence alignments")
plt.ylabel("Number of trials")
plt.xlabel("Alignment score")
plt.show()
print(" and bar plotting.\n")

mean = 0
for key in scoring_distribution.keys():
    mean += scoring_distribution[key] * key
mean /= iterations

standard_deviation = 0
for key in scoring_distribution.keys():
    for dummy_iterations in range(key):
        standard_deviation += (scoring_distribution[key] - mean) ** 2
standard_deviation = math.sqrt(standard_deviation / iterations)

z_score = (hf_local_alignment[0] - mean) / standard_deviation

print("Question 5: Mean =", mean, "\nStandard deviation =",
      standard_deviation, "\nZ-score =", z_score, "\n")

print("Question 6: Based on the Z-score found in question 5, it is clear",
      "that the local alignment for human and fruitfly proteins is not",
      "caused by chance. In comparison, winning the lottery is more likely",
      "to happen than the proteins being a product of chance.\n")

print("Question 7:\ndash score for edit distance = 0\ndiagonal score for edit",
      "distance = 2\noff diagonal score for edit distance = 1\n")

print("Question 8:", end=" ")
spelling_list = read_words("Data//Word_list.txt")
humble_spelling = check_spelling("humble", 1, spelling_list)
firefly_spelling = check_spelling("firefly", 2, spelling_list)
print("Words with edit distance = 1 from Humble:", humble_spelling,
      "\n\nWords with edit distance = 2 from Firefly:", firefly_spelling)
