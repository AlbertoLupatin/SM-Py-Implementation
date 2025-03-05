import argparse
import pandas as pd
import numpy as np

# Check if the input strings are DNA sequences
def CheckDNA(s):
    dna = {"A", "C", "T", "G"}
    if not set(s).issubset(dna):
        raise ValueError(f"Sequence {s} is not a valid DNA sequence")
    return s

# Compute the scoring matrix
def scoringMatrix(s1, s2, match, mismatch, gap):
    sub = pd.DataFrame(np.zeros((len(s1), len(s2))), index=s1, columns=s2, dtype=int)
    path = pd.DataFrame("None", index=s1, columns=s2, dtype=object)

    for r in range(1, len(s1)):
        for c in range(1, len(s2)):
            mat_mis = sub.iloc[r-1, c-1] + (match if s1[r] == s2[c] else mismatch)
            gap_up = sub.iloc[r-1, c] + gap
            gap_left = sub.iloc[r, c-1] + gap

            sub.iloc[r, c] = max(mat_mis, gap_up, gap_left, 0)
            
            if sub.iloc[r, c] == mat_mis:
                path.iloc[r, c] = ("match" if s1[r] == s2[c] else "mismatch")
            elif sub.iloc[r, c] == gap_up:
                path.iloc[r, c] = "up"
            elif sub.iloc[r, c] == gap_left:
                path.iloc[r, c] = "left"
            else:
                path.iloc[r, c] = "None"
    
    path.iloc[0,:] = "None"
    path.iloc[:,0] = "None"

    return sub, path

# Perform traceback to get the best alignment
def perform_traceback(sub, path, start_r, start_c):
    l1, l2 = [], []
    curr_r, curr_c = start_r, start_c
    curr = path.iloc[curr_r, curr_c]

    while curr != "None":
        if curr == "match" or curr == "mismatch":
            l1.append(path.index[curr_r])
            l2.append(path.columns[curr_c])
            curr_r, curr_c = curr_r-1, curr_c-1

        elif curr == "up":
            l1.append(path.index[curr_r])
            l2.append("-")
            curr_r -= 1

        elif curr == "left":
            l1.append("-")
            l2.append(path.columns[curr_c])
            curr_c -= 1

        curr = path.iloc[curr_r, curr_c]

    l1.reverse()
    l2.reverse()
    return "".join(l1), "".join(l2), sub.iloc[start_r, start_c]

# Get the best alignment
def Traceback(s1, s2, sub, path):
    r, c = np.unravel_index(np.argmax(sub.to_numpy()), sub.shape)
    res1, res2, score = perform_traceback(sub, path, r, c)
    return res1, res2, score

# Run the Smith-Waterman algorithm
def SmithWaterman(s1, s2, matScore, misScore, gapPen):
    s1 = CheckDNA(s1.upper())
    s2 = CheckDNA(s2.upper())
    
    s1, s2 = [""] + list(s1), [""] + list(s2)
    sub, path = scoringMatrix(s1, s2, matScore, misScore, gapPen)
    
    res1, res2, score = Traceback(s1, s2, sub, path)
    print(f"Best Alignment\n{'-'*40}\nSequence 1: {res1}\nSequence 2: {res2}\nScore: {score}")

# Argument parsing
parser = argparse.ArgumentParser(
    description=(
        "------------------------------------------------------------\n"
        "Smith-Waterman Algorithm - Local Sequence Alignment\n"
        "Author: Alberto Lupatin\n"
        "Algorithmics for Bioinformatics Exam, University of Trento\n"
        "------------------------------------------------------------\n"),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("s1", type=str, help="First DNA sequence")
parser.add_argument("s2", type=str, help="Second DNA sequence")
parser.add_argument("--match_score", type=int, default=3, help="Score for matching bases. Default: 3")
parser.add_argument("--mismatch_score", type=int, default=-3, help="Penalty for mismatching bases. Default: -3")
parser.add_argument("--gap_pen", type=int, default=-2, help="Gap penalty. Default: -2")

args = parser.parse_args()
SmithWaterman(args.s1, args.s2, matScore=args.match_score, misScore=args.mismatch_score, gapPen=args.gap_pen)