#Results

Across microbial datasets, SAD predictions from the maximum entropy theory of ecology (METE) generally failed to explain more than 60% of variation in abundance among species (Figure 1; Table 1). This is a relatively poor level of explanatory power given that METE commonly explains 90% or more of variation among macroscopic plants and animals.

Previous tests of METE on communities of macroscopic plants and animals **(are much more successful)** produced much greater success, often explaining 90% or greater variation in abundance (Harte et al. 2008, 2009, White et al. 2012, Xiao et al. 2014).  
**S3:** As expected, due to its relatively even form, the broken-stick model (effectively the geometric distribution) performed considerably worse than METE and generally produced negative r-square values.   
**S4:** Negative values were possible because the relationship is not fitted, i.e., estimating variation around a line with a slope of 1.0 and intercept of zero (White et al. 2012, Locey and White 2013, Xiao et al. 2014).  
**S5:** While the log-series (METE) characterizes the form of the SAD better than the broken-stick, microbial SADs are still characterized by disparities in abundance that METE fails to capture.  
**S6:** Both METE and the broken-stick under-predict the abundance of the most abundant species and over-predict the abundance of the rarest species.

##P2
**S1:** We found that the success of METE and the broken-stick were influenced by the two primary state-variables (*N* and *S*) and the primary constraint of average abundance (*N*/*S*) (Table 1). Across each dataset (EMP, HMP, MG-RAST) increasing *N* led to decreasing fits of each model while 




Table 1. 

| Dataset      | Model | Variable |  r<sup>2</sup>  | p-value |
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
|   MGRAST 95% |  BS   |     N    | -0.302 | 0.141
|   MGRAST 95% |  METE |     N    | -0.158 | 0.828
|   MGRAST 95% |  BS   |     S    | 0.0234 | 0.828
|   MGRAST 95% |  METE |     S    | 0.140 | 0.192
|   MGRAST 95% |  BS   |     N/S  | -0.862 | 3.75*10<sup>-27</sup> 
|   MGRAST 95% |  METE |     N/S  |-0.734 | 4.12*10<sup>-16</sup> 
|   MGRAST 97% |  BS   |     N    | -0.0782 | 0.480
|   MGRAST 97% |  METE |     N    | 0.226 | 0.0389
|   MGRAST 97% |  BS   |     S    | 0.169 | 0.125
|   MGRAST 97% |  METE |     S    | 0.353 | 0.00101
|   MGRAST 97% |  BS   |     N/S  | -0.642 | 4.69*10<sup>-11</sup>
|   MGRAST 97% |  METE |     N/S  | -0.244 | 0.0255
|   MGRAST 99% |  BS   |     N    | -0.312 | 0.00265
|   MGRAST 99% |  METE |     N    | -0.172 | 0.109
|   MGRAST 99% |  BS   |     S    | 0.0150 | 0.890
|   MGRAST 99% |  METE |     S    | 0.132 | 0.221
|   MGRAST 99% |  BS   |     N/S  | -0.868 | 7.99*10<sup>-28</sup>
|   MGRAST 99% |  METE |     N/S  | -0.737 | 2.71*10<sup>-16</sup>
