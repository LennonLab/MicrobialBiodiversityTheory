# Discussion
Within and among communities of macroscopic organisms, METE often explains 90% or more of observed variation in abundance among species. 
Here, we showed that while METE performs better than an alternative MaxEnt prediction (i.e., Broken-stick) it often fails to explain the majority of variation within and among communities of bacteria and archaea. 
These results are primarily due to the tendency of both models to under-predict dominance the abundance of the most abundant species.

We also showed that METE's success is heavily influenced by one of its primary state variables (*N*). 
As a result, increasing *N* causes METE as well as the Broken Stick to fail more severely. 
Importantly, these conditions also characterize numerical differences between microbial and macrobial SAD datasets. 
That is, *N* for microbial datasets often represents tens of thousands to millions of processed rRNA reads  while *N* for macrobial SADs typically ranges from a few hundred to a few thousand individual organisms. Consequently, METE might fail for microbes because it can be expected to fail with increasing *N*.

The failure of the Broken-stick model and METE could have been anticipated. It has been shown that as *N* increases, the evenness of the SAD can be expected to decrease as a result of numerical constraints (Locey and White 2013). 
In the same way, as average abundance (*N*/*S*) increases, the evenness of the SAD can be expected to naturally decrease. 
In both cases, constraints on the form of the SAD imposed by *N* and *N*/*S* lead to increasingly uneven SADs that outstrip the relatively even form predicted by the Broken-stick (i.e. the geometric distribution) as well as the relatively uneven form predicted by METE (i.e. the log-series distribution). 
Still, it remains to be seen whether the inability of METE to predict microbial SADs is entirely driven by numerical constraints.

Our study suggests that highly uneven SADs are driven by factors leading to high *N*. 
However, uneven microbial SADs could also be driven by factors suggested to explain the microbial rare biosphere. 
For example, widespread dispersal and the ability of microbes to persist in suboptimal environments may allow many small populations of dormant or slow-growing organisms to persist (Reid and Buckley 2011).  Additionally, microorganisms seem to have unparalleled capacities to partition limited resources (Muscarella et al. 2014). Consequently, the failure of the Broken-stick and METE may owe as much to the statistical influence of *N* as to the ecological mechanisms that cause differences in abundances among specific species.

Our study reveals that ecology may lack an SAD model that accurately predicts differences in abundance among microbes. 
More generally ecology may lack an appropriate model to predict abundances when *N* scales beyond a few tens of thousands. 
At the least, ecology currently lacks a MaxEnt based SAD model that doesn't fail with increasing *N*.
To this end, we suggest an SAD model with naturally greater unevenness than the log-series distriubtion, i.e., the Zipf-distribution (Zipf 1949). The Zipf-distribution is based on a power-law and characterized by a single free parameter:

$$N(k) \propto \exp (-bk)/k^{\gamma}$$

and can be derived as a prediction of MaxEnt (Baek et al. 2011, Visser 2013). Likewise, the Zipf-distribution has been previously identified as the SAD model that best fits microbial SADs (Gans 2005). 
However, this particular finding has received little attention.

## Conclusion
Constraint-based biodiversity theory provides a first-principle framework for predicting biodiversity patterns based solely on empirical inputs. 
While MaxEnt frameworks infer no ecological processes they also avoid a common source of error in the study of commonness and rarity, i.e., inferring the importance of processes that have gone unverified and unmeasured (e.g., speciation, niche differentiation, coexistence). 
Yet, it is clear from our study that the maximum entropy theory of ecology (METE) will fail for communities of very large *N*, such as microbiomes where *N* easily exceeds many million. 
Consequently, while the microbial "rare biosphere" appears to be an exceptional phenomenon, we cannot yet say whether the cause is due more to the biology driving rarity independent of *N* or to the biology that allows microbes to attain such high degrees of *N*. 
But what we can say is that microbial ecology and perhaps macroecology lacks an SAD theory to accurately predict distributions of abundance among large ecological communities, microbial communities, or both.
