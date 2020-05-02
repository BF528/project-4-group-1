# Contents of the folder
1. Sample Metadata.txt : Metadata of the sample analysed in the project
2. whitelist.qsub: shell script to create whitelist of barcodes to be used during alignment
3. alevin.qsub: script to run salmon alevin to get cell-by-gene matrix for each library
4. mergematrix.py: script to merge the 3 cell-by-gene matrices into one UMI counts matrix for further processing
