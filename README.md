# MicroMETE
Project repository for data and Python code associated with the testing of maximum entropy theory in microbial ecology. 


## Re-running the Analyses/ Generating the Figures

Files in the data folder ending with ".zip" will have to be unzipped before the code can be run.

The code accepts the following arguments from the user. 

**Flags**

**-f:** If the data is already in the file path, a figure will be generated using this argument. The argument specifies which figure is generated. 

| Argument |          Figure         |
|:--------:|:-----------------------:|
|     1    |         Figure 1        |
|     2    |         Figure 2        |
|    S1    |  Supplementary figure 1 |
|    S2    |  Supplementary figure 2 |
|    S3    |  Supplementary figure 3 |
|    S4    |  Supplementary figure 4 |
|    S5    |  Supplementary figure 5 |
| S6       | Supplementary figure 6  |
| S7       | Supplementary figure 7  |
| S8       | Supplementary figure 8  |
| S9       | Supplementary figure 9  |
| S10      | Supplementary figure 10 |

 
**-a:**  Argument for rerunning the analysis used to generate the observed vs. predicted file and NSR2 file for a given figure. Inputs are 'Yes' or 'No.' Default is 'No.'


**-r:** Run all analyses and generate all figures used in the main body and supplement of the paper.  



**Ex:** To generate the analysis and figure for supplemental figure 7, the script would be run with the following arguments.  

	python generate_figs_zipf.py -a yes -f S7

## Dependencies

Python version 2.7.10-2 is used. 

The following Python modules are used in this analysis.

+ numpy 1.10.1-py27_0

+ matplotlib 1.4.2-np19py27_0

+ scipy 0.16.0

+ [macroeco_distributions](https://github.com/weecology/macroecotools)

+ [macroecotools](https://github.com/weecology/macroecotools)

+ [mete](https://github.com/weecology/METE)

## The MIT License (MIT)

Copyright (c) 2015  William R. Shoemaker, Kenneth J. Locey, Jay T. Lennon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

## Attributes

To write the code in this repository we used MIT liscensed code from the GitHub repositories [METE](https://github.com/weecology/macroecotools) and [macroecotools](https://github.com/weecology/macroecotools) on 9/20/2015. 
