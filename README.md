# MicrobialBiodiversityTheory
Project repository for data and Python code associated with the following preprint:

Shoemaker WR, Locey KJ, Lennon JT. (2016) A unifying ecological theory of microbial biodiversity. PeerJ Preprints 4:e1450v3 https://doi.org/10.7287/peerj.preprints.1450v3

## Re-running the analyses and re-generating the figures

Files in the data folder ending with ".zip" will have to be unzipped before the code can be run.

To Figure2 and supplementary  figures 1 and 2 you will need to request data, as the files are too large to store on GitHub. The data will need to be unzipped and in the main folder in a sub-folder named "data."

The code accepts the following arguments from the user.

**Flags**

**-a** or **--analysis:** Runs the analysis required for the figure.

**-f:**  The figure you want to generate. Indicate what figure you want to generate after the flag using the table below.

| Argument |          Figure         |
|:--------:|:-----------------------:|
|     1    |         Figure 1        |
|     2    |         Figure 2        |
|    3    |  Figure 3  |
|    4    | Figure 4 |
|    S1    |  Supplementary figure 1 |
|    S2    |  Supplementary figure 2 |


### Order of operations

If you want to regenerate the data (warning, this is very computationally intensive and will take several days to complete) start with step 1, otherwise run step 2 for the figures.

1) Run the following command.

	python runAnalysis.py -a

2) Run this command

	python runAnalysis.py -f

  and indicate what figure  you want after "-f"


## Dependencies

Python version 2.7.10-2 is used.

The following Python modules/versions are used in this analysis.

+ numpy 1.10.1

+ matplotlib 1.4.2

+ Pandas 0.16.2

+ scipy 0.16.0

+ setuptools 18.4

+ statsmodels 0.6.1

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
