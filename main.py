'''
ALGORITHM FOR OPTIMAL PAIRWISE SEQUENCE ALIGNMENT
USING SEQUENCE EXAMPLE FROM INTRO TO BIOINFORMATICS ATHUR M.LESK PAGE 139 EXAMPLE 4.5.

SCORING SCHEME:
    MATCH = 0
    MISMATCH = 20
    INSERTION / DELETION = 25

USES A MINIMUM FUNCTION FOR WEIGHTED DISTANCE.
'''


from useful_funcs import score

# example sequences
init_seq2 = "ATG"
init_seq1 = "GGAATGG"

gap_penalty = 25

n = len(init_seq1)
m = len(init_seq2)

dp = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

for i in range(1, n + 1):
    dp[i][0] = i * gap_penalty  # across i the gaps

for j in range(1, m + 1):
    dp[0][j] = j * gap_penalty  # down j and adding gaps

# for row in dp:
#     print(row)

# fill out the matrix

for i in range(1, n + 1):
    for j in range(1, m + 1):
        dp[i][j] = min(
            dp[i - 1][j - 1] + score(init_seq1[i - 1], init_seq2[j - 1]),
            dp[i - 1][j] + gap_penalty,
            dp[i][j - 1] + gap_penalty
        )
for row in dp:
    print(row)


# perform traceback to find the optimal sequence

i, j = n, m

seq_1 = ''
seq_2 = ''

while i > 0 or j > 0:
    # get current position in the matrix
    curr = dp[i][j]

    # match / mismatching base pair. Go diagonal.
    if i > 0 and j > 0 and curr == (dp[i - 1][j - 1] + score(init_seq1[i - 1], init_seq2[j - 1])):
        # this is a match so add to both strings.
        seq_1 = init_seq1[i - 1] + seq_1
        seq_2 = init_seq2[j - 1] + seq_2
        i -= 1
        j -= 1

    # goes horizontal - therefore gap in sequence 1 (seq one going along left of matrix)
    elif j > 0 and curr == dp[i][j - 1] + gap_penalty:
        seq_1 = '-' + seq_1
        seq_2 = init_seq2[j - 1] + seq_2
        j -= 1

    # goes up vertically - therefore gap in sequence 2 (which is going along the top of the matrix)
    elif i > 0 and curr == (dp[i - 1][j] + gap_penalty):
        seq_1 = init_seq1[i - 1] + seq_1
        seq_2 = '-' + seq_2
        i -= 1

print('Aligned Sequences:')
print(seq_1)
print(seq_2)


