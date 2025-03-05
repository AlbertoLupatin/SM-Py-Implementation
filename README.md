# Smith-Waterman Algorithm for Local Sequence Alignment  

## 📌 Project Overview  
This project implements the **Smith-Waterman algorithm**, a dynamic programming approach for **local sequence alignment**.  
It identifies the most similar subsequences between two given DNA, RNA, or protein sequences, using **affine gap penalties** and a substitution matrix.  

This project was developed as part of the *Algorithmics for Bioinformatics* course at the **University of Trento**.  

## 🚀 Features  
- Implements the **Smith-Waterman algorithm** for local alignment.  
- Supports **custom scoring matrices** (e.g., PAM, BLOSUM).  
- Uses **affine gap penalties** for improved accuracy.  
- Outputs the **optimal local alignment** with a traceback mechanism.  
- Efficient implementation with **dynamic programming**.  

## 🛠️ Installation & Requirements  
This project requires **Python 3.x** and the following dependencies:  

```bash
pip install numpy
```

## ▶️ Usage

Run the script with:
```
python smith_waterman.py -s1 SEQUENCE1 -s2 SEQUENCE2 -m MATRIX -g GAP_PENALTY
```

**Parameters**
	•	-s1 and -s2: Input sequences (DNA, RNA, or Protein).
	•	-m: Scoring matrix (PAM, BLOSUM, or a custom file).
	•	-g: Gap penalty (negative value).

## 📚 References
	•	Smith, T. F., & Waterman, M. S. (1981). Identification of common molecular subsequences. Journal of Molecular Biology, 147(1), 195–197.
	•	Bioinformatics Algorithms - University of Trento Lecture Notes.
