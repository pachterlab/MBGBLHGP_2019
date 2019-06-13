# Code for reproducing Figure 1, Supplementary Figures 1, and Supplementary Note Figure 2 from the paper "Modular and efficient pre-processing of single-cell RNA-seq"

After pre-processing your data and making a corrected and sorted bus file, then `bustools text -o output.txt output.correct.sort.bus` will make a text bus file. Then run make_bug.py to make a bug.txt file. All scripts process this bug file.
