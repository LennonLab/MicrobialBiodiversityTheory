#Results
##P1
**S1:** Across microbial datasets, SAD predictions from the maximum entropy theory of ecology (METE) generally failed to explain more than 60% of variation in abundance among species and often produced negative r-squared values, which have no straigtforward interpretation (Figure 1; Table 1).  
**S2:** Previous tests of METE on communities of macroscopic plants and animals **(are much more successful)** produced much greater success, often explaining 90% or greater variation in abundance (Harte et al. 2008, 2009, White et al. 2012, Xiao et al. 2014).  
**S3:** As expected, due to its relatively even form, the broken-stick model (effectively the geometric distribution) performed considerably worse than METE and generally produced negative r-square values.   
**S4:** Negative values were possible because the relationship is not fitted, i.e., estimating variation around a line with a slope of 1.0 and intercept of zero (White et al. 2012, Locey and White 2013, Xiao et al. 2014).  
**S5:** While the log-series (METE) characterizes the form of the SAD better than the broken-stick, microbial SADs are still characterized by disparities in abundance that METE fails to capture.  
**S6:** Both METE and the broken-stick under-predict the abundance of the most abundant species and over-predict the abundance of the rarest species.

##P2
**S1:** We found that the success of METE and the broken-stick were influenced by the two primary state-variables (*N* and *S*) and the primary constraint of average abundance (*N*/*S*) (Table 1). Across each dataset (EMP, HMP, MG-RAST) increasing *N* led to decreasing fits of each model while 




Table 1. 

| Dataset      | Model | Variable |  r<sup>2</sup>  | p-value |
|:------------:|:-----:|:--------:|:-----:|:-------:|
|   HMP        |   BS  |     N    |-0.386 |   1.15*E<sup>-159</sup>   |
|   HMP        |  METE |     N    |-0.191 |   2.01*E<sup>-38</sup>   |
|   HMP        |   BS  |     S    | 0.276 |   2.82*E<sup>-79</sup>      |
|   HMP        |  METE |     S    | 0.314 |   1.44*E<sup>-103</sup>       |
|   HMP        |   BS  |    N/S   |-0.626 |   0.0   |
|   HMP        |  METE |    N/S   |-0.453 |   1.87*E<sup>-226</sup>       |
|   EMP closed |   BS  |     N    |-0.354 |   0.0   |
|   EMP closed |  METE |     N    |-0.0824| 2.02*E<sup>-23</sup> |
|   EMP closed |   BS  |     S    | 0.264 |  4.89*E<sup>-231</sup>       |
|   EMP closed |  METE |     S    | 0.287 |1.32*E<sup>-274</sup>        |
|   EMP closed |   BS  |    N/S   |-0.695 |   0.0   |
|   EMP closed |  METE |    N/S   |-0.377 |   0.0   |
|   EMP open   |  BS   |    N     |-0.349 |   0.0   |
|   EMP open   |  METE |    N     |-0.205 |   6.28*E<sup>-140</sup>      |
|   EMP open   |  BS   |    S     | 0.0731| 5.00*E<sup>-19</sup>         |
|   EMP open   |  METE |    S     | 0.103 | 1.57*E<sup>-36</sup>        |
|   EMP open   |  BS   |    N/S   |-0.763 | 0.0      |
|   EMP open   |  METE |    N/S   |-0.544 | 0.0      |
