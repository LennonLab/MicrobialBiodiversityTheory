# Introduction

Understanding patterns of abundance, distribution, and diversity is a central goal of ecology (Brown 1995). 
Among the most universal of these patterns is the observation that few species in most ecological communities are relatively abundant, while most are relatively rare. 
Explaining this has been a major focus of community ecology and ecological theory for several decades (e.g., Motomura 1932, Fisher et al. 1943, Preston 1948, MacArthur 1957, Whittaker 1972, Sugihara 1980, Tokeshi 1990, Hubbell 2001, Tilman 2004). 
But, while many models predict uneven species abundance distributions (SAD) as the result of resource partitioning, dispersal limitation and demographic stochasticity, and competition and coexistence, the most successful models have been statistical (e.g., Fisher et al. 1943, Preston 1948).

A new paradigm in predicting the SAD takes a statistical constraint-based approach (Harte et al. 2008, 2009, Harte 2011, Pueyo et al. 2006, White et al. 2012, Locey and White 2013, Xiao et al. 2014). 
Recognizing that the form of the SAD is constrained by total abundance (*N*) and the number of species (*S*), several constraint-based models predict the most likely form of the SAD given *N* and *S* as the only empirical inputs (e.g., Harte et al. 2008, 2009, Harte 2011, Pueyo et al. 2006, Haegeman and Etienne 2010). 
Most of these constraint-based approaches use a maximum-likelihood framwork that invokes the principle of maximum entropy, i.e., where the most likely form of the SAD is that having the most ways of occurring given *N*, *S*, and/or average abundance (*N*/*S*).

Though simple in concept, MaxEnt SAD predictions require assumptions such as whether species and individuals are distinguishable (Pueyo et al. 2006, Harte 2011). Consequently, different predictions based on *N* and *S* are possible under a MaxEnt framework (Haegeman and Etienne 2010). So far, however, the maximum entropy theory of ecology (METE) has been the most accurate MaxEnt framework for predicting the SAD and in linking it to other primary patterns of biodiversity (Harte 2011, Xiao et al. 2014).

METE has been widely successful in predicting the SAD among communities of macroscopic organisms (White et al. 2012). METE often explains more than 90% of observed variation in abundance among species of plant and animal, based on the largest available compilations of ecological community data (White et al. 2012). Yet, despite its success in predicting SADs and attempts to unify other biodiversity patterns through the SAD, METE has yet been tested among the most abundant and taxonomically diverse organisms on Earth, i.e., bacteria and archaea.

Within natural and host-associated ecosystems, most microbial taxa account for the minority of relative abundance. 
This seemingly universal pattern of microbial commonness and rarity is known as the microbial "rare biosphere". 
Though often described, e.g., as the percent of species with less than 0.1%  of total abundance, the rare biosphere pattern is rarely predicted with general models of biodiversity (e.g. neutral theory, niche-theory, METE). 
Yet, the rare biosphere pattern appears to reiterate the universal nature of uneven SADs, the most powerful predictions of which are often produced by METE.  

Here, we test the SAD prediction of METE using the largest compilation of microbial community data ever assembled from publicly available sources. 
These data include 20,216 sites of bacterial and archaeal communities from the Earth Microbiome Project, the Human Microbiome Project, and datasets from Argonne National Laboratory's metagenomic server MG-RAST. 
As an additional MaxEnt prediction and comparison to METE, we use the prediction for the Broken Stick model, which is also based on *N* and *S* (Haegeman and Etienne 2010).
In demonstrating previously undocumented failure of METE to predict SADs, we suggest a more appropriate MaxEnt model for SADs characterized by large *N*, i.e., the MaxEnt form of the Zipf-distribution.
