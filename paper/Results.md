#Results

SAD predictions from the maximum entropy theory of ecology (METE) generally explained 0 to less than 55% of variation in HMP and EMP data (Figure 1; Table 1).
This is a poor level of explanatory power given that METE commonly explains 90% or more of variation among macroscopic plants and animals (Baldridge 2015). METE, however, performed considerably better for MG-RAST datasets, often explaining 81% of variation among microbes (Figure 1; Table 1). 

Differences in the performance of METE among datasets are well-explained by differences in *N*.
We found that MG-RAST data were characterized by smaller values of *N* (and the best performance of METE), while EMP open-reference data were characterized by the highest values of *N* and the worst performance of METE.
Across all datasets, the success of METE and the Broken-stick were influenced by *N*, where increasing *N* led to decreasing performance of each model (Table 2).

The percent sequence similarity cutoff used to cluster 16S rRNA reads into operational taxonomic units had no effect on the explanatory of METE and the Broken-Stick, even though this should influence the value of *S*. However, we did find that the amount of the variation explained by the geometric distribution increased with *S*, an expected result that has been previously predicted (Wilson 1993).


Table 1.

| Dataset | Model | $$\overline{r^{2}_{m}}$$ | $$\sigma_{\bar{r^{2}}}$$ | $$N$$ | $$S$$ |
|:--------:|:-----:|:-------:|:------:|:----:|:-----:|
|  HMP     |  BS | -0.543  | 0.0170  | 5050 | 78 |
|  HMP     |  METE   | 0.520   | 0.00846 ||
|EMP closed|  BS | -0.434  | 0.00851 | 44779 | 1189 |
|EMP closed|  METE   | 0.562   | 0.00377 ||
|EMP open  |  BS | -0.881  | 0.0101  | 88751 | 7247 |
|EMP open  |  METE   | 0.0619  | 0.00526 ||
|MGRAST 95%|  BS | 0.551   | 0.0242  | 1200 | 247|
|MGRAST 95%|  METE   | 0.816   | 0.0113  ||
|MGRAST 97%|  BS | 0.571   | 0.0184  | 929 | 210|
|MGRAST 97%|  METE   | 0.816   | 0.0103  ||
|MGRAST 99%|  BS | 0.542   | 0.0244  | 1148 | 235 |
|MGRAST 99%|  METE   | 0.811   | 0.0113  ||




Table 2. 

| Dataset      | Model | Variable |  $r$  | p-value |
|:------------:|:-----:|:--------:|:-----:|:-------:|
|   HMP        |   BS  |     N    |-0.386 |   1.15*10<sup>-159</sup>   |
|   HMP        |  METE |     N    |-0.191 |   2.01*10<sup>-38</sup>   |
|   HMP        |   BS  |     S    | 0.276 |   2.82*10<sup>-79</sup>      |
|   HMP        |  METE |     S    | 0.314 |   1.44*10<sup>-103</sup>       |
|   HMP        |   BS  |    N/S   |-0.626 |   0.0   |
|   HMP        |  METE |    N/S   |-0.453 |   1.87*10<sup>-226</sup>       |
|   EMP closed |   BS  |     N    |-0.354 |   0.0   |
|   EMP closed |  METE |     N    |-0.0824| 2.02*10<sup>-23</sup> |
|   EMP closed |   BS  |     S    | 0.264 |  4.89*10<sup>-231</sup>       |
|   EMP closed |  METE |     S    | 0.287 |1.32*10<sup>-274</sup>        |
|   EMP closed |   BS  |    N/S   |-0.695 |   0.0   |
|   EMP closed |  METE |    N/S   |-0.377 |   0.0   |
|   EMP open   |  BS   |    N     |-0.349 |   0.0   |
|   EMP open   |  METE |    N     |-0.205 |   6.28*10<sup>-140</sup>      |
|   EMP open   |  BS   |    S     | 0.0731| 5.00*10<sup>-19</sup>         |
|   EMP open   |  METE |    S     | 0.103 | 1.57*10<sup>-36</sup>        |
|   EMP open   |  BS   |    N/S   |-0.763 | 0.0      |
|   EMP open   |  METE |    N/S   |-0.544 | 0.0      |
|   MGRAST 95% |  BS   |     N    | -0.302 | 0.141 |
|   MGRAST 95% |  METE |     N    | -0.158 | 0.828 |
|   MGRAST 95% |  BS   |     S    | 0.0234 | 0.828 |
|   MGRAST 95% |  METE |     S    | 0.140 | 0.192 |
|   MGRAST 95% |  BS   |     N/S  | -0.862 | 3.75*10<sup>-27</sup> | 
|   MGRAST 95% |  METE |     N/S  |-0.734 | 4.12*10<sup>-16</sup> | 
|   MGRAST 97% |  BS   |     N    | -0.0782 | 0.480 |
|   MGRAST 97% |  METE |     N    | 0.226 | 0.0389 |
|   MGRAST 97% |  BS   |     S    | 0.169 | 0.125 |
|   MGRAST 97% |  METE |     S    | 0.353 | 0.00101 |
|   MGRAST 97% |  BS   |     N/S  | -0.642 | 4.69*10<sup>-11</sup> |
|   MGRAST 97% |  METE |     N/S  | -0.244 | 0.0255 |
|   MGRAST 99% |  BS   |     N    | -0.312 | 0.00265 |
|   MGRAST 99% |  METE |     N    | -0.172 | 0.109 |
|   MGRAST 99% |  BS   |     S    | 0.0150 | 0.890 |
|   MGRAST 99% |  METE |     S    | 0.132 | 0.221 |
|   MGRAST 99% |  BS   |     N/S  | -0.868 | 7.99*10<sup>-28</sup> |
|   MGRAST 99% |  METE |     N/S  | -0.737 | 2.71*10<sup>-16</sup> |
