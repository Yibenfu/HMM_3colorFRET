# HMM_3colorFRET
Hidden Markov model (HMM) with the constraint EM algorithm for the analysis of 3color FRET data. 
The algorithm is coded in MATLAB. It can read in the original data, analyze it with 6 states (5 on DNA and 1 in solution), and output the analyzed results.
To run the code, please follow the below steps:

1. download all the files in one folder.
2. put the original data in the same folder.
3. open maincode.m file by MATLAB.
4. At the top of the maincode.m, in the file-read commend line, change the file's name according to your data file.
5. run the code. The fitted data will be saved in a folder with the folder name the same as the original data file's name.

Note: one may need to adjust initial approximations of FRET values, variance values, and transition matrix. 
