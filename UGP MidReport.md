# Mid-Term Report: Meta-GWAS Analysis 

#### Dhruv Gupta
#### **Date:** March 2026

This is a Mid-Term Report detailing the tasks given and the progress on them / approach used.

---
# Task 1: *Generating a simulated GWAS summary stat and plot Manhattan, Q-Q and volcano plots*

The first task was to  simulate a GWAS dataset using the package BioGWAS and then generate summary statistics.
  
### 1.1 What is GWAS?

- Genome Wide Association Studies (GWAS) are studies that deal with scanning the entire human genome to identify certain genetic variants called Single Nucleotide Polymorphisms or SNPs associated with specific phenotypes (could look like traits like height, weight, skin colour etc. or disease risk like cholestrol levels or response to drugs etc.)
- For example, when studying creatinine levels (a kidney function marker), GWAS helps identify which genetic variants make some people naturally have higher or lower creatinine.

#### Notebook 1: `bioGWAS_simulated.ipynb` - Basic GWAS Simulation and Summary Statistics

### 1.2 What is BioGWAS?

[Link to bioGWAS paper here](https://www.mdpi.com/2079-7737/13/1/10)
What the BioGWAS framework does:
- simulates genotype data matrix (N*M) under some frequency of allele assumptions
- defines a subset of causal gene variants
- generates phenotypes under specific genetic architecture
- performs SNP-wise association

**What was done:
- Simulated 10,000 SNPs with random chromosomal positions
- Generated genotype data (0, 1, 2 copies of effect allele) for 2,000 individuals
- Created a simulated quantitative trait with 20 causal SNPs having real effects
- Performed association testing using linear regression
- Generated standard GWAS summary stats (QQ, Manhattan and Volcano Plots)

Note: **Linkage Disequilibrium was set to zero for this set of simulations, ie, the SNPs are all assumed to be independent in this notebook.**

The model to simulate data is:
$$ y = X \beta + \epsilon $$
- X: Genotype matrix
- $\beta$: SNP effects (sparse vector, mostly zero since most SNPs dont effect)
- $\epsilon$: noise term for simulation

Assumptions:
- SNPs all have only two alleles
- hence, genotypes are only encoded as 0, 1, or 2 depending on effected allele count
- some independence assumption amongst SNPs -> **ASSUMING NO LINKAGE DISEQUILIBRIUM**

![[Pasted image 20260313225340.png]]
Fig. 1: Manhattan Plot for simulated BioGWAS data with zero LD 
  
![[Pasted image 20260313225434.png]]
Fig 2: Volcano plot for simulated BioGWAS data with zero LD 


#### Notebook 2: `bioGWAS_with_LD.ipynb` - GWAS Data Simulation with Non-Zero LD and Summary Statistics

### 1.3 Linkage Disequilibrium

**Linkage Disequilibrium (LD)** is the non-random association of alleles at different loci in a population. 
- SNPs that are physically close on a chromosome tend to be inherited together
- This creates **correlation** between nearby genetic variants
- LD decays with physical distance due to recombination
- Measured by $r^2$ score usually

Note: LD matters for GWAS because a detected association signal might not always be due to a causal variant, it could be due to a correlated one in linkage diseq. with it. Hence, ignoring LD will lead to unrealistic simulations.

### 1.4 Simulation Model

The phenotype model is still:
$$ y = X \cdot \beta + \epsilon $$
**difference:** The genotype matrix $X$ is now generated with correlation structure:
- SNPs close together have higher correlation
- LD decays exponentially with distance
- $\text{Cor}(\text{SNP}_i, \text{SNP}_j) = \rho^{|i-j|}$ where $\rho$ is the LD decay parameter

![[Pasted image 20260313225850.png]]
Fig. 3: Linkage Disequilibrium Heatmap depicting the correlation between various SNPs

![[Pasted image 20260313225927.png]]
Fig 4: Manhattan Plot for BioGWAS simulated data with Non-Zero LD

![[Pasted image 20260313230006.png]]
Fig 5: Volcano Plot for BioGWAS simulated data with Non-Zero LD

---
# Task 2: *To perform MetaGWAS analysis on Simulated Datasets*

#### Notebook 3: `bioGWAS_metaGWAS.ipynb` - Meta-GWAS done using 4 individual GWAS studies 
### 2.1 Why Meta-Analysis?

Individual GWAS studies often have limited statistical power due to small sample sizes. **Meta-GWAS** combines results from multiple studies to:
- Increase statistical power to detect true genetic associations
- Validate findings across different populations
- Identify associations that individual studies might miss

### 2.2 Meta-GWAS

**Challenges with combining multiple studies:**
- Studies may have different sample sizes
- Different ancestry backgrounds (allele frequencies vary)
- Different genotyping platforms
- Heterogeneous effect sizes

**Two Main Approaches:**
- **Fixed-Effects (FE)**: Assumes all studies estimate the same true effect
- **Random-Effects (RE)**: Allows true effects to vary across studies

#### Fixed-Effects Model 
(weighted using Inverse-Variance)
Combines effect estimates weighted by their precision:
$$\beta_{\text{meta}} = \frac{\sum_{i=1}^{k} w_i \beta_i}{\sum_{i=1}^{k} w_i}$$

where $w_i = \frac{1}{SE_i^2}$ (inverse variance weight)
Standard error:
$$SE_{\text{meta}} = \sqrt{\frac{1}{\sum_{i=1}^{k} w_i}}$$

#### Heterogeneity Assessment
**Cochran's Q statistic** tests if effect sizes are homogeneous:
$$Q = \sum_{i=1}^{k} w_i (\beta_i - \beta_{\text{meta}})^2$$

**I² statistic** quantifies proportion of variance due to heterogeneity:
$$I^2 = \frac{Q - (k-1)}{Q} \times 100\%$$

- $I^2 < 25\%$: Low heterogeneity
- $I^2 = 25-75\%$: Moderate heterogeneity  
- $I^2 > 75\%$: High heterogeneity

We simulated 4 different independent GWAS studies using BioGWAS: (shown in table below):

1. **Different sample sizes** (500, 1000, 1500, 2000 participants)
2. **Different ancestry-like characteristics** (varying MAF distributions)
3. **Same true causal variants** (but different power to detect them)
4. **Linkage disequilibrium** (realistic correlation structure)
5. **Some heterogeneity** (effect sizes vary slightly across studies)

![[Pasted image 20260313230301.png]]
Fig 6: Table depicting the simulated 4 GWAS studies to be used for Meta-GWAS.

![[Pasted image 20260313230417.png]]
Fig 7: Manhattan plots of all individual GWAS studies followed by MetaGWAS study combining them. It is clear that increase in number of SNPs leads to increase of chance of identifying causal SNPs.


---
# Task 3: *To perform MetaGWAS analysis on ACTUAL Datasets* - In Progress.

This is a task that is currently in progress. So far, I have switched to using built-in frameworks like METAL package to do this due to better reproducibility. The task was for us to validate our results with results from the lab (ie, validate the plots etc.) 

- Creatinine was chosen as the trait of choice to perform MetaGWAS using the METAL framework.
- The following plots were plotted (shown below):
	- Volcano Plot
	- Manhattan Plot
	- QQ Plot

These plots are pasted below and are currently in the process of verification. 

![[Pasted image 20260313231149.png]]
Fig 8: Volcano Plot for Creatinine, showing significant SNPs for the real dataset.


![[Pasted image 20260313231255.png]]
Fig. 9: Manhattan Plot for Creatinine, showing the address of most significant SNP.


![[Pasted image 20260313231324.png]]
Fig 10: QQ Plot for Creatine.


