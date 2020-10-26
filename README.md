# DNA sequence alignment analysis

## Matrices and functions
For this project , two types of matrices will be used: alignment matrices and scoring matrices. The first two functions that you will implement compute a common class of scoring matrices and compute the alignment matrix for two provided sequences, respectively. The first function builds a scoring matrix as a dictionary of dictionaries.

+ `build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score)`: Takes as input a set of characters `alphabet` and three scores `diag_score`, `off_diag_score`, and `dash_score`. The function returns a dictionary of dictionaries whose entries are indexed by pairs of characters in `alphabet` plus `’-’`. The score for any entry indexed by one or more dashes is `dash_score`. The score for the remaining diagonal entries is `diag_score`. Finally, the score for the remaining off-diagonal entries is `off_diag_score`.

One final note for `build_scoring_matrix` is that, although an alignment with two matching dashes is not allowed, the scoring matrix should still include an entry for two dashes (which will never be used). The second function computes an alignment matrix using the method `ComputeGlobalAlignmentScores`. The function computes either a global alignment matrix or a local alignment matrix depending on the value of `global_flag`.

+ `compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag)`: Takes as input two sequences `seq_x` and `seq_y` whose elements share a common alphabet with the scoring matrix `scoring_matrix`. The function computes and returns the alignment matrix for `seq_x` and `seq_y`. 

If `global_flag` is `True` the algorithm will use:

![CGAS algorithm](https://github.com/rokobo/DNA-Sequence-Analysis/blob/main/Data/Compute%20Global%20Alignment%20Scores.jpg?raw=True)

If `global_flag` is `False`, each entry is computed using the method described in above, but with the modification: Whenever Algorithm `ComputeGlobalAlignmentScores` computes a value to assign to `S[i,j]`, if the computed value is negative, the algorithm instead assigns 0 to `S[i,j]`.

## Alignment functions
For the second part, the alignment matrix returned by `compute_alignment_matrix` will be used to compute global and local alignments of two sequences `seq_x` and `seq_y`. The first function will implement the method ComputeAlignment:

![CA algorithm](https://github.com/rokobo/DNA-Sequence-Analysis/blob/main/Data/Compute%20Alignment.jpg?raw=True)

+ `compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix)`: Takes as input two sequences (`seq_x` and `seq_y`) whose elements share a common alphabet with the scoring matrix `scoring_matrix`. This function computes a global alignment of `seq_x` and `seq_y` using the global alignment matrix `alignment_matrix`.The function returns a tuple of the form `(score, align_x, align_y)` where `score` is the score of the global alignment `align_x` and `align_y`. Note that `align_x` and `align_y` should have the same length and may include the padding character `’-’`.

This second function will compute an optimal local alignment starting at the maximum entry of the local alignment matrix and working backwards to zero:

+ `compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix)`: Takes as input two sequences (`seq_x` and `seq_y`) whose elements share a common alphabet with the scoring matrix `scoring_matrix`. This function computes a local alignment of `seq_x` and `seq_y` using the local alignment matrix `alignment_matrix`.The function returns a tuple of the form `(score, align_x, align_y)` where `score` is the score of the optimal local alignment `align_x` and `align_y`. Note that `align_x` and `align_y` should have the same length and may include the padding character `’-’`.

## Protein comparison
In following questions, we compute the similarity between the human and fruit fly versions of the eyeless protein and see if we can identify the PAX domain.

### Question 1 
Load the files `HumanEyelessProtein` and `FruitflyEyelessProtein`. These files contain the amino acid sequences that form the eyeless proteins in the human and fruit fly genomes, respectively. Then load the scoring matrix `PAM50` for sequences of amino acids. This scoring matrix is defined over the alphabet `{A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, Z, X, -}` which represents all possible amino acids and gaps (the "dashes" in the alignment). Next, compute the local alignments of the sequences of `HumanEyelessProtein` and `FruitflyEyelessProtein` using the `PAM50` scoring matrix in order to find the score and local alignments for these two sequences. 

### Question 2 
To continue our analysis, we next consider the similarity of the two sequences in the local alignment computed in Question 1 to a third sequence. The file `ConsensusPAXDomain` contains a "consensus" sequence of the PAX domain; that is, the sequence of amino acids in the PAX domain in any organism. In this problem, we will compare each of the two sequences of the local alignment computed in Question 1 to this consensus sequence to determine whether they correspond to the PAX domain.

Load the file ConsensusPAXDomain. For each of the two sequences of the local alignment computed in Question 1, do the following:

+ Delete any dashes `’-’` present in the sequence.

+ Compute the global alignment of this dash-less sequence with the `ConsensusPAXDomain` sequence.

+ Compare corresponding elements of these two globally-aligned sequences (local vs. consensus) and compute the percentage of elements in these two sequences that agree.

To reiterate, you will compute the global alignments of local human vs. consensus PAX domain as well as local fruitfly vs. consensus PAX domain. Your answer should be two percentages: one for each global alignment.

### Question 3 
Examine your answers to Questions 1 and 2. Is it likely that the level of similarity exhibited by the answers could have been due to chance? In particular, if you were comparing two random sequences of amino acids of length similar to that of `HumanEyelessProtein` and `FruitflyEyelessProtein`, would the level of agreement in these answers be likely? Include a short justification.

## Hypothesis testing
One weakness of our approach in Question 3 was that we assumed that the probability of any particular amino acid appearing at a particular location in a protein was equal. In the next two questions, we will consider a more mathematical approach to answering Question 3 that avoids this assumption. In particular, we will take an approach known as statistical hypothesis testing to determine whether the local alignments computed in Question 1 are statistically significant.

### Question 4 
Write a function `generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials)` that takes as input two sequences `seq_x` and `seq_y`, a scoring matrix `scoring_matrix`, and a number of trials `num_trials`. This function should return a dictionary `scoring_distribution` that represents an un-normalized distribution generated by performing the following process `num_trials` times:

+ Generate a random permutation `rand_y` of the sequence `seq_y` using `random.shuffle()`.

+ Compute the maximum value `score` for the local alignment of `seq_x` and `rand_y` using the score matrix `scoring_matrix`.

+ Increment the entry `score` in the dictionary `scoring_distribution` by one.

Use the function `generate_null_distribution` to create a distribution with 1000 trials using the protein sequences `HumanEyelessProtein` and `FruitflyEyelessProtein` (using the PAM50 scoring matrix). Next, create a bar plot of the normalized version of this distribution. The horizontal axis should be the scores and the vertical axis should be the fraction of total trials corresponding to each score. 

### Question 5 
Given the distribution computed in Question 4, we can do some very basic statistical analysis of this distribution to understand how likely the local alignment score from Question 1 is. To this end, we first compute the mean μ and the standard deviation σ of this distribution via:

<img src="https://render.githubusercontent.com/render/math?math=\displaystyle \mu = \frac{1}{n}\sum_i s_i">

<img src="https://render.githubusercontent.com/render/math?math=\displaystyle \sigma = \sqrt{\frac{1}{n}\sum_i(s_i - \mu)^2}">

where the values <img src="https://render.githubusercontent.com/render/math?math=s_i"> are the scores returned by the `n` trials. If `s` is the score of the local alignment for the human eyeless protein and the fruitfly eyeless protein, the z-score `z` for this alignment is:

<img src="https://render.githubusercontent.com/render/math?math=z = \frac{s-\mu}{\sigma}">

The z-score helps quantify the likelihood of the score ss being a product of chance. Small z-scores indicate a greater likelihood that the local alignment score was due to chance while larger scores indicate a lower likelihood that the local alignment score was due to chance. Calculate:

+ The mean and standard deviation for the distribution that you computed in Question 4. 

+ The z-score for the local alignment for the human eyeless protein vs. the fruitfly eyeless protein based on these values.
### Question 6 
For bell-shaped distributions such as the normal distribution, the likelihood that an observation will fall within three multiples of the standard deviation for such distributions is very high. Based on your answers to Questions 4 and 5, is the score resulting from the local alignment of the `HumanEyelessProtein` and the `FruitflyEyelessProtein` due to chance? As a concrete question, which is more likely: the similarity between the human eyeless protein and the fruitfly eyeless protein being due to chance or winning the jackpot in an extremely large lottery? Provide a short explanation.

## Spelling correction
In other applications, measuring the dissimilarity of two sequences is also useful. Given two strings, the edit distance corresponds to the minimum number of single character insertions, deletions, and substitutions that are needed to transform one string into another. In particular, if `x` and `y` are strings and aa and bb are characters, these edit operations have the form:

+ Insert - Replace the string `x+y` by the string `x+a+y`.

+ Delete - Replace the string `x+a+y` by the string `x+y`.

+ Substitute - Replace the string `x+a+y` by the string `x+b+y`.

### Question 7 
Similarity between pairs of sequences and edit distances between pairs of strings are related. In particular, the edit distance for two strings `x` and `y` can be expressed in terms of the lengths of the two strings and their corresponding similarity score as follows:

<img src="https://render.githubusercontent.com/render/math?math=\text{Similarity = }\mid x\mid + \mid y \mid -\text{score}(x, y)">

Where `score(x, y)` is the score returned by the global alignment of these two strings using a very simple scoring matrix that can be computed using `build_scoring_matrix`. Also, determine the values for `diag_score`, `off_diag_score`, and `dash_score` such that the score from the resulting global alignment yields the edit distance when substituted into the formula above. Indicating which values corresponds to which parameters.

### Question 8 
In practice, edit distance is a useful tool in applications such as spelling correction and plagiarism detection where determining whether two strings are similar/dissimilar is important. Implement a simple spelling correction function that uses edit distance to determine whether a given string is the misspelling of a word. To begin, load the word_list of 79339 words. Then, write a function `check_spelling(checked_word, dist, word_list)` that iterates through word_list and returns the set of all words that are within edit distance `dist` of the string `checked_word`. Use your function `check_spelling` to compute the set of words within an edit distance of one from the string "humble" and the set of words within an edit distance of two from the string "firefly".