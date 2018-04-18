# MicrobialBiodiversityTheory
Project repository for data and Python code associated with the following preprint:

Shoemaker WR, Locey KJ, Lennon JT. (2016) A unifying ecological theory of microbial biodiversity. PeerJ Preprints 4:e1450v3 https://doi.org/10.7287/peerj.preprints.1450v3

## Re-running the analyses and re-generating the figures

The following files and folders need to be uncompressed before you run the analsis:

`~/MicrobialBiodiversityTheory/data/EMPclosed-Data/EMPclosed-SADs.txt.zip`

`~/MicrobialBiodiversityTheory/data/EMPopen-Data/EMPopen-SADs.txt.zip`

`~/MicrobialBiodiversityTheory/data/HMP-Data/hmp1.v35.hq.otu.counts.bz2`

`~/MicrobialBiodiversityTheory/data/ObsPred.zip`


To Figure2 and supplementary figures 1 and 2 you will need to request data, as the files are too large to store on GitHub.

The code accepts the following arguments from the user.

**Flags**

**`-a`** or **`--analysis`:** Runs the analysis required for the figure.

**`-f`:**  The figure or table you want to generate. Indicate what figure/ table you want to generate after the flag using the table below.

| Argument |          Figure/ Table         |
|:--------:|:-----------------------:|
|     F1    |         Figure 1        |
|     F2    |         Figure 2        |
|    F3    |  Figure 3  |
|    F4    | Figure 4 |
|    FS1    |  Supplementary figure 1 |
|    FS2    |  Supplementary figure 2 |
|    FS3    |  Supplementary figure 3 |
|    FS4    |  Supplementary figure 4 |
|    FS5    |  Supplementary figure 5 |
|    FS6    |  Supplementary figure 6 |
|     T1    |         Table 1        |
|     T2    |         Table 2        |
|    TS1    |  Supplementary table 1 |
|    TS2    |  Supplementary table 2 |
|    TS3    |  Supplementary table 3 |
|    TS4    |  Supplementary table 4 |
|    TS5    |  Supplementary table 5 |
|    TS6    |  Supplementary table 6 |

### Order of operations

If you want to regenerate the data (warning, this is very computationally intensive and will take several days to complete) start with step 1, otherwise run step 2 for the figures.

1) Run the following command.

	`python runAnalysis.py -a`

2) Run this command

	`python runAnalysis.py -f`

  and indicate what figure  you want after "`-f`"


## Dependencies

`Python 2.7.10-2` is used.

The following Python modules/versions are used in this analysis.

+ `numpy 1.10.1`

+ `matplotlib 1.4.2`

+ `Pandas 0.16.2`

+ `scipy 0.16.0`

+ `setuptools 18.4`

+ `statsmodels 0.6.1`

+ [`macroecotools 0.2`](https://github.com/weecology/macroecotools)

+ [`mete 0.1`](https://github.com/weecology/METE)


## A note of caution

Since we wrote the code for this analysis (~2015) the [Weecology](http://www.weecology.org/) group has updated [`macroecotools`](https://github.com/weecology/macroecotools). Newer versions of macroecotools are incompatible with this code. In addition, as described in the manuscript, we ran into issues with fitting the lognormal using maximum likelihood estimation on microbial communities with a large number of individuals. If you want to just reproduce our results, then using [`macroecotools 0.2`](https://github.com/weecology/macroecotools) is fine. **However, if you want to fit the lognormal to your data, we strongly recommend that you work with the most recent version of [`macroecotools`](https://github.com/weecology/macroecotools)**. 


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

To write the code in this repository we used MIT liscensed code from the GitHub repositories [mete 0.1](https://github.com/weecology/METE) and [macroecotools 0.2](https://github.com/weecology/macroecotools) on 9/20/2015.
