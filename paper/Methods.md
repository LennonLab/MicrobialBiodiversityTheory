# Methods

## Data
For the analysis we used bacterial and archaeal community sequence data from 15,535 sites.  14,962 of these sites were from the Earth Microbiome Project (EMP) (Gilbert et al., 2014) obtained on 22 August, 2014. Sample processing and sequencing of the V4 region of the 16s ribosomal RNA gene are standardized by the EMP and all are publicly available at www.microbio.me/emp. The EMP data consist of open and closed reference datasets, which are defined in the QIIME tutorial (http://qiime.org/tutorials/otu_picking.html).  
QIIME defines closed-reference as a classification scheme where any reads that do not hit a sequence in a reference collection are excluded from analysis. In contrast, open-reference refers to a scheme where reads that do not hit a reference collection are subsequently clustered de novo and represent unique but unclassified taxonomic units. Our main results are based on closed-reference data, due to the greater accuracy of the approach and because unclassified sequences were excluded from other microbial datasets (below).	
We also used 4,303 sites from the Data Analysis and Coordination Center (DACC) for the National Institutes of Health (NIH) Common Fund supported Human Microbiome Project (HMP). These data consisted of samples taken from 15 to 18 locations (including the skin, gut, vagina, and oral cavity) on each of 300 healthy individuals. 
In each sample the V3-V5 region of the 16S rRNA gene was sequenced and analyzed using the mothur pipeline (Turnbaugh, et al., 2007). We excluded sites from pilot phases of the HMP as well as time-series data; see http://hmpdacc.org/micro_analysis/microbiome_analyses.php. for details on HMP sequencing and sampling protocols.

We also included 1,319 non-experimental sequencing projects consisting of processed 16s rRNA amplicon reads from the Argonne National Laboratory metagenomics server MG-RAST (Meyer, et al., 2008).  

Represented in this compilation were samples from arctic aquatic systems (130 sites; MG-RAST id: mgp138), hydrothermal vents (123 sites; MG-RAST id: mgp327) (Flores et al., 2011), freshwater lakes in China (187 sites; MG-RAST id: mgp2758) (Wang, et al., 2014), arctic soils (44 sites; MG-RAST id: mgp69) (Chu et al., 2010), temperate soils (84 sites; MG-RAST id: mgp68) (Fierer et al., 2012), bovine fecal samples (16 sites; MG-RAST id: mgp14132), human gut microbiome samples not part of the HMP project (529 sites; MG-RAST id: mgp401) (Yatsunenko, et al., 2012), a global-scale dataset of indoor fungal systems (128 sites) (Amend et al., 2010), and freshwater, marine, and intertidal river sediments (34 sites; MG-RAST id: mgp1829). 


A common convention in lieu of traditional species classificaiton for microbial community sequence data is to cluster 16s rRNA amplicon reads into Operational Taxonomic Units (OTUs) based on a 97% cutoff for sequence similarity. Locey and White showed that the percent cutoff of sequence similarity does not change the shape of the SAD (Locey & While, 2013). However, how the percent cutoff affects the fit of SAD models to emperical data is rarely been tested in the literature (Dumbrell et al., 2010; Woodcock et al., 2007). The use of MG-RAST allowed us to choose common parameter values for percent sequence similarity (i.e. 97% for species-level) and taxa assignment including a maximum e-value (probability of observing an equal or better match in a database of a given size) of 10-5, a minimum alignment length of 50 base pairs, and minimum percent sequence similarities of 95, 97, and 99% to the closest reference sequence in MG-RAST’s M5 rRNA database (Flores et al., 2011; Wang, et al., 2014; Chu et al., 2010; Fierer et al., 2012; Yatsunenko, et al., 2012; Amend et al., 2010). Quantifying dominance, evenness, rarity, and richness. We calculated or estimated aspects of diversity (dominance, evenness, rarity, richness) for each site in our data compilation. All analyses can be reproduced or modified for further exploration by using code, data, and following directions provided here: https://github.com/LennonLab/MicroMETE.  

## MaxEnt predictions of the SAD
### METE
The maximum entropy theory of ecology (METE) (Harte et al. 2008, 2009, Harte 2011) is based on two empirical inputs: species richness (*S*) and total abundance (*N*). These, along with an inferred rate of community-level metabolism (*E*), form the state variables of METE. Four constraints are produced from these state variables. These are the average number of individuals per species (*N*/*S*), the average per species metabolic flux (*E*/*S*), and the constraints that no species has more than *N* individuals or a greater total metabolic rate than *E*. *E* is later integrated out of the SAD prediction.  

The prediction of METE is based on a joint conditional probability distribution that describes the distribution of individuals (*n*) over species and of metabolism (*ε*) over individuals within a species (Harte et al. 2008, Harte 2011). Entropy of the distribution is then maximized according to the method of Lagrangian multipliers (Jaynes 2003, Harte 2011). The SAD is then derived by integrating out energy and dropping terms that are vanishingly small. This process then yields the log-series SAD (Fisher et al. 1943). The log-series distribution is among the oldest and most successful SAD models but has generally lacked a convincing first-principle explanation from either an ecological or statistical perspective. In this case, METE predicts the shape of which is dependent only on the values of *S* and *N*:

$$\Phi\left ( n\mid S_{0},N_{0} \right ) = \frac{1}{log(\beta ^{-1})}\frac{e^{-\beta n}}{n}
$$

where $$\beta$$ is defined by the equation 

$$\frac{N_{0}}{S_{0}}=\frac{\sum_{n=1}^{N_{0}}e^{-\beta n}}{\sum_{n=1}^{N_{0}}e^{-\beta n}/n}
$$

### Broken-stick 
 While some other MaxEnt models produce similar, if not, identical (Pueyo et al. 2007, Dewar and Porté 2008, Frank 2011) predictions for the SAD, MaxEnt models based on different assumptions can yield very different predictions (Haegeman and Etienne 2010). One example is the simultaneous discrete Broken-stick model of MacArthur (1960), which as pointed out by Haegeman and Etienne (2010) is simply the geometric distribution with mean *N*/*S*. Unlike the log-series, the broken-stick model predicts a relatively even distribution which is often a poor fit to empirical SADs (Hubbell 2001). The broken-stick gives equal weight to all ordered configurations of *S* species whose abundances sum to *N*, the equation for which for the $$r^{th}$$ rarest species being:

$$\frac{N}{S}\sum_{i=1}^{r}\frac{1}{S-i+1}$$

With $$r$$ being the abundance of the $$r^{th}$$
rarest species. 

## Testing MaxEnt predictions
Both METE (which predicts a log-series distribution) and the Broken-stick (i.e., the geometric distribution) produce predictions for the rank-abundance form of the SAD. This form of the SAD is simply a vector of species abundances ranked from greatest to least. Both predictions yield the same value of *S* that is given as the empirical input. This means that the observed and predicted SADs can be directly compared using regression analyses to reveals the percent variation explained by each model (METE, Broken-stick). We generated the predicted forms of the SAD using the source code of White et al. (2012) (https://github.com/weecology/white-etal-2012-ecology) and the macroecotools repository (https://github.com/weecology/macroecotools), which contains functions for fitting maximum-likelihood forms of species abundance models in addition to functions for other macroecological analyses. Using that source code, we calculated the modified coefficient of determination ($$r_{m}^{2}$$) around the 1-to-1 line (as per White et al. 2012, Locey and White 2013, Xiao et al. 2014).

$$r_{m}^{2} = 1 - \frac{sum((obs - pred)^{2})}{sum((obs-(\overline{obs}))^{2})}$$

Negative values were possible because the relationship is not fitted, i.e., estimating variation around a line with a slope of 1.0 and intercept of zero (White et al. 2012, Locey and White 2013, Xiao et al. 2014).