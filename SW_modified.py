import argparse
import pandas as pd
import numpy as np

# Check if the input strings are DNA sequences
def CheckDNA(s):
    dna = {"A", "C", "T", "G"}
    if not set(s).issubset(dna):
        raise ValueError(f"Sequence {s} is not a valid DNA sequence")
    return s

def scoringMatrix(s1, s2, match, mismatch, gap):
    s1, s2 = [""] + list(s1), [""] + list(s2)

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

# Filter the sequences by the optional parameters
def checkParameters(s1, s2, curr_sc, score_req, curr_gap, gap_req, min_occ):
    if curr_sc < score_req:
        return False

    if curr_gap < gap_req:
        return False
    
    for i in range(len(s1) - min_occ + 1):
        if s1[i : i+min_occ] == s2[i : i+min_occ]:
            return True

    return False

def getMaxScore(sub):
    max_score = sub.to_numpy().max()
    (r, c) = np.unravel_index(np.argmax(sub.to_numpy()), sub.shape)

    return max_score, (r, c)

# Perform the actual traceback with the given values
def perform_traceback(sub, path, start_r, start_c):
    l1, l2 = [], []
    curr_r, curr_c = start_r, start_c
    curr = path.iloc[curr_r, curr_c]

    matCount, misCount, gapCount = 0, 0, 0

    while curr != "None":
        if curr == "match" or curr == "mismatch":
            l1.append(path.index[curr_r])
            l2.append(path.columns[curr_c])
            if curr == "match":
                matCount += 1
            else:
                misCount += 1
            curr_r, curr_c = curr_r-1, curr_c-1

        elif curr == "up":
            l1.append(path.index[curr_r])
            l2.append("-")
            gapCount += 1
            curr_r -= 1

        elif curr == "left":
            l1.append("-")
            l2.append(path.columns[curr_c])
            gapCount += 1
            curr_c -= 1

        curr = path.iloc[curr_r, curr_c]

    l1.reverse()
    l2.reverse()
    return "".join(l1), "".join(l2), sub.iloc[start_r, start_c], gapCount

# Exploits the perform_traceback function to return the traceback outputs
def Traceback(s1, s2, sub, path, score_req, gap_req, ret_max, min_occ):
    if ret_max == 0:
        bestScore, (r, c) = getMaxScore(sub)
        res1, res2, _, _ = perform_traceback(sub, path, r, c)
        
        return res1, res2, bestScore

    elif ret_max == 1:
        all_res = {}

        bestScore, _ = getMaxScore(sub)

        score_req = 0.01 * score_req * bestScore

        for r in range(2, len(s1)):
            for c in range(2, len(s2)):
                if path.iloc[r, c] == "None":
                    continue
                
                l1, l2, curr_sc, gapCount = perform_traceback(sub, path, r, c)
    
                consecutive_check = False
                for i in range(len(s1) - min_occ + 1):
                    if s1[i : i+min_occ] == s2[i : i+min_occ]:
                        consecutive_check = True

                if curr_sc >= score_req and gapCount >= gap_req and consecutive_check:
                    all_res[(l1, l2)] = (curr_sc, len(l1), gapCount)

        return dict(sorted(all_res.items(), key=lambda item: item[1][1], reverse=True))

    else:
        raise ValueError("ret_max option must be 0 (best score) or 1 (all alignments)")

# Print the filtered sequences
def printSequences(sequences):
    print("Alignments found:")
    print("-"*40)
    for (seq1, seq2), (score, length, gap) in sequences.items():
        print(f"Sequence 1: {seq1}\nSequence 2: {seq2}\nScore: {score}\nLength: {length}\tGaps: {gap}\n{'-'*40}")

# Exploits all the function above to compute the Smith-Waterman algorithm
def SmithWaterman(s1, s2, matScore, misScore, gapPen, ret_max, score_req, gap_req, min_occ):
    s1 = CheckDNA(s1.upper())
    s2 = CheckDNA(s2.upper())

    sub, path = scoringMatrix(s1, s2, matScore, misScore, gapPen)

    if ret_max == 0:
        res1, res2, score = Traceback(s1, s2, sub, path, score_req, gap_req, ret_max, min_occ)
        print(f"Best Sequence\n{'-'*40}\nSequence 1: {res1}\nSequence 2: {res2}\nScore: {score}")

    else:
        res = Traceback(s1, s2, sub, path, score_req, gap_req, ret_max, min_occ)
        printSequences(res)


# Decription of the program
parser = argparse.ArgumentParser(
    description=(
        "------------------------------------------------------------\n"
        "Smith-Waterman Algorithm - Local Sequence Alignment\n"
        "Author: Alberto Lupatin\n"
        "Algorithmics for Bioinformatics Exam, University of Trento\n"
        "------------------------------------------------------------\n"),
    formatter_class=argparse.RawTextHelpFormatter
                
)

# Input values
parser.add_argument("s1", type = str, help = "First DNA sequence")
parser.add_argument("s2", type = str, help = "Second DNA sequence")

parser.add_argument("--match_score", type = int, help = "Score for matching bases (default: 3)", default=3)
parser.add_argument("--mismatch_score", type = int, help = "Penalty for mismatching base (default: -3)", default=-3)
parser.add_argument("--gap_pen", type = int, help = "Gap penalty (default: -2)", default=-2)
parser.add_argument("--ret_max", type = int, help = "0: Returns only the best alignment\n1 (default): Returns all the alignments with the chosen filters", default=1)

parser.add_argument("--score_req", type = int, help = "Minimum alignment score required - set as percentage of the maximum score (default: 60)", default=60)
parser.add_argument("--gap_req", type = int, help = "Minimum number of gaps allowed (default: 1)", default=1)
parser.add_argument("--min_occ", type = int, help="Number of consecutive matches involving the same nucleotide (default: 3)", default=3)


args = parser.parse_args()
SmithWaterman(args.s1, args.s2, matScore=args.match_score, misScore=args.mismatch_score, gapPen=args.gap_pen, 
              ret_max=args.ret_max, score_req=args.score_req, gap_req=args.gap_req, min_occ=args.min_occ)
