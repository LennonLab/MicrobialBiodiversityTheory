# Methods

## Data


## MaxEnt predictions of the SAD
### METE
**P1**  
**S1:** The maximum entropy theory of ecology (METE) (Harte et al. 2008, 2009, Harte 2011) is based on two empirical inputs: species richness (*S*) and total abundance (*N*).   
**S2:** These, along with an inferred rate of community-level metabolism (*E*), form the state variables of METE.  
**S3:** Four constraints are produced from these state variables.  
**S4:** These are the average number of individuals per species (*N*/*S*), the average per species metabolic flux (*E*/*S*), and the constraints that no species has more than *N* individuals or a greater total metabolic rate than *E*.  
**S5:** *E* is later integrated out of the SAD prediction.  

**P2**  
**S1:** The prediction of METE is based on a joint conditional probability distribution that describes the distribution of individuals (*n*) over species and of metabolism (*ε*) over individuals within a species (Harte et al. 2008, Harte 2011). 
**S2:** Entropy of the distribution is then maximized according to the method of Lagrangian multipliers (Jaynes 2003, Harte 2011).  
**S3:** The SAD is then derived by integrating out energy and dropping terms that are vanishingly small. This process then yields the log-series SAD (Fisher et al. 1943).  
**S4:** The log-series distribution is among the oldest and most successful SAD models but has generally lacked a convincing first-principle explanation from either an ecological or statistical perspective.  
**S4:** In this case, METE predicts the shape of which is dependent only on the values of *S* and *N*:

###insert equation 1

where β is defined by the equation 

###insert equation 2


### Broken-stick 
**S1:** While some other MaxEnt models produce similar, if not, identical (Pueyo et al. 2007, Dewar and Porté 2008, Frank 2011) predictions for the SAD, MaxEnt models based on different assumptions can yield very different predictions (Haegeman and Etienne 2010).   
**S2:** One example is the geometric distribution, which as pointed out by Haegeman and Etienne (2010) is also the simultaneous discrete Broken-stick model of MacArthur (1960).  
**S3:** Unlike the log-series, the broken-stick model predict a relatively even distribution that is often a poor fit to empirical SADs (ref).  
**S4:** The broken-stick gives equal weight to all ordered configurations of *S* species whose abundances sum to *N*, the equation for which is:

###insert equation 3

## Testing MaxEnt predictions
**S1:** Both METE (log-series) and the Broken-stick (geometric distribution) produce predictions for the rank-abundance form of the SAD.  
**S2:** That is, vectors of species abundances ranked from greatest to least.  
**S3:** Both predictions yield the same value of *S* that is given as the empirical input.  
**S4:** This means that the observed and predicted SADs can be directly compared using regression analyses to reveals the percent variation explained by each model (METE, Broken-stick).  
**S5:** We found the predicted form of the SAD for METE and the broken-stick, for each site in our data, then calculated a modified coefficient of determination (r-square) around the 1-to-l line (as per White et al. 2012, Locey and White 2013, Xiao et al. 2014).
 