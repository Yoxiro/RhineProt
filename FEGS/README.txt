The algorithm of FEGS is implemented by Matlab language, and we have prepared a test sample.

---Usage---
Before executing FEGS, users should set path firstly as follows:
$ addpath PATH-TO-FEGS\ .
$ Put the data to be analyzed in the current folder.

Execute FEGS as follows:
$ FEGS(fasta_data); 

---Input and Output files---
Input: fasta_data is the input file name of fasta type containing multiple protein sequences. 
Output: the matrix of n rows and m columns, where n is the number of protein sequences, and m is the number of extracted features.

---Example---
To run the test data, users can use the following command:
$ FEGS ('data')