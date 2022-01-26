Peptide intensity input file 
=============================

What is expected?
------------------

If you upload your own data, we expect a text file which contains at least the protein names, peptide sequences, modifications per peptide and intensities. If multiple samples were run in the same experiment, each sample should have its own intensity column. 

Example
-----------

| Protein 	|       Peptide    	|     Modification 	|    Intensity_sampleA 	|   Intensity_sampleB 	|
|:-------:	|:-------------:	|:-------------:	|:-------:	|:-------:	|
|  P16043 	| RHVDAIFTTNYRK 	| [5] Oxidation 	|  11247  	|  10378  	|

Options to read in file
-----------------------

- **header**: is there a header (line with column names, generally first line) present in your data file? Check box if present.
- **Number of lines to skip**: are there text lines present above the header that should not be read in? If yes, enter how many such lines are present.
- **File separator**: which file separator separates the different columns in your file? Check the right one.
- **Intensity columns identifier**: string that is present in all columns that contain the intensity values. In the example above this would be "Intensity".
- **Protein column**: The number of the column that contains the protein ID. In the example above this would be 1.
- **Sequence column**: The number of the column that contains the peptide sequence. In the example above this would be 2.
- **Modifications column**: The number of the column that contains the peptide modification. In the example above this would be 3.

Once you have filled out all options correctly, click on the button "Read in your own data". 