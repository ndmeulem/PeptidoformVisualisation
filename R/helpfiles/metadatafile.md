Metadata input file
======================

What is expected?
------------------------

We expect a file that contains the names of the intensity columns in one column and in the other columns identifying information. These columns should contain at least enough information so that each row represents one unique ID (without taking the intensity column into account).


Example 
-----------

| Intensity   columns 	|   species 	  |   treatment     	|   sample 	  |
|:-------------------:	|:-------:	|:---------:	|:------:	|
| 210112_hESCs_DDA_P2 	|  human  	|     S     	|    1   	|
| 210114_hESCs_DDA_P3 	|  human  	|     S     	|    2   	|
| 210114_hESCs_DDA_P4 	|  human  	|     S     	|    3   	|
| 210115_hESCs_DDA_P6 	|  human  	|     S     	|    4   	|
| 210115_hESCs_DDA_P8 	|  human  	|     NS     	|    5   	|
| 210114_hESCs_DDA_P5 	|  human  	|     NS     	|    6   	|
| 210115_hESCs_DDA_P7 	|  human  	|     NS     	|    7   	|
| 210118_hESCs_DDA_P1 	|  human  	|     NS     	|    8   	|
| 210118_hESCs_DDA_P9 	|  human  	|     NS     	|    9   	|
| 210112_hESCs_DDA_N1 	|  mouse  	|     S     	|   10   	|
| 210112_hESCs_DDA_N2 	|  mouse  	|     S     	|   11   	|
| 210114_hESCs_DDA_N3 	|  mouse  	|     S     	|   12   	|
| 210114_hESCs_DDA_N4 	|  mouse  	|     S     	|   13   	|
| 210115_hESCs_DDA_N5 	|  mouse  	|     NS     	|   14   	|
| 210115_hESCs_DDA_N6 	|  mouse  	|     NS     	|   15   	|
| 210115_hESCs_DDA_N7 	|  mouse  	|     NS     	|   16   	|
| 210115_hESCs_DDA_N8 	|  mouse  	|     NS     	|   17   	|
| 210118_hESCs_DDA_N9 	|  mouse  	|     NS     	|   18   	|

Options 
---------

- **File separator**: which file separator separates the different columns in your file? Check the right one.