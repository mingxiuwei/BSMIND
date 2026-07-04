## BSMIND 

Current version: 0.1.1 (2026-07-02)  


The R package BSMIND is developed for the integrative analysis of spatial multi-omics data to identify tissue spatial domains. 
It jointly models transcriptomic and proteomic data to capture consensus spatial structures, 
which can effectively handle complete and incomplete domain configurations across two omics datasets. 
With omics expression matrices and spatial coordinates as basic inputs, users can conveniently complete multi-omics integrative clustering, 
quantitative performance evaluation of clustering results, and follow-up biological enrichment analysis. 


## Prerequisites and Installation

1. R version >= 4.1.3

2. CRAN packages: MCMCpack, truncnorm, mvtnorm

3. Install the package `BSMIND`

```R
devtools::install_github("mingxiuwei/BSMIND")
```


## Datasets information

The data description is given in the following table.

|               Dataset           | Spot/Cell number |  Gene number |  ADT number  |                        Download links                        |
| :-----------------------------: | :--------------: | :----------: | :----------: | :------------------------------------------------------------------------: |
|      Human lymph node A1        |       3,484      |    18,085    |      31      | Raw data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263617     |
|  MMTV-PyMT mouse breast cancer  |       1,978      |    18,932    |      32      | Raw data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi                   |


## Example Code


### 1. simulation 1

The following code shows an example (the first simulation study in the manuscript) that runs the function "BSMIND" in our package.

Import the required R package.


```R
library(BSMIND)
```

Load the example data stored in this package. The example data includes:

data_1: The simulated transcriptomic data.
data_2: The simulated proteomic data.
coord: The simulated coordinates (2 columns).
true_label: The true consensus labels in this simulation study.


```R
sim1_1_data <- load_sim1_1_data()
data_1 <- sim1_1_data$data1
data_2 <- sim1_1_data$data2
coord <- sim1_1_data$coord
true_label <- sim1_1_data$true_label
```

Run "BSMIND" function to obtained consensus labels. The meaning of each argument in the function is listed below.

* data_1: First n*p1 data matrix (spots x genes), gene expression modality.
* data_2: Second n*p2 data matrix (spots x features), secondary feature modality.
* coord: Data frame of spot coordinates, 2 columns: x, y coordinates.
* platform: Spatial sequencing platform; "ST" (square grid) or "Visium" (hex grid).
* K_fixed: Fixed number of regions for all latent layers. Default is 3.
* a_mu: Mean of normal prior for gene expression cluster means mu_gk. Default is 0.
* b_mu: Standard deviation of normal prior for mu_gk. Default is 1.
* IGkappa: Shape parameter of inverse-gamma prior for gene expression sigma_g. Default is 2.
* IGtau: Scale parameter of inverse-gamma prior for sigma_g. Default is 10.
* a_mu2: Mean of normal prior for second modality cluster means mu_pk. Default is 2.
* b_mu2: Standard deviation of normal prior for mu_pk. Default is 1.
* IGkappa2: Shape parameter of inverse-gamma prior for second modality sigma_p. Default is 2.
* IGtau2: Scale parameter of inverse-gamma prior for sigma_p. Default is 10.
* a_beta: Prior mean for spatial Potts interaction parameter beta. Default is 0.7.
* tau_beta: Prior standard deviation for beta. Default is 1.
* dpAlpha: Dirichlet concentration parameter for cluster proportion prior pi_vec. Default is 1.
* minPsi1: Lower bound of uniform prior for gene expression contribution weight psi1. Default is 0.8.
* maxPsi1: Upper bound of uniform prior for gene expression contribution weight psi1. Default is 0.95.
* minPsi2: Lower bound of uniform prior for secondary feature modality contribution weight psi2. Default is 0.8.
* maxPsi2: Upper bound of uniform prior for secondary feature modality contribution weight psi2. Default is 0.95.
* tau0: Scale parameter for Dirichlet proposal distribution of pi_vec. Default is 2.
* tau1: Proposal standard deviation for Potts parameter beta. Default is 0.05.


Obtain consensus clustering labels (C_new) for Simulation 1.
```R
set.seed(51686)
result <- BSMIND(
  data_1 = data_1,
  data_2 = data_2,
  coord = coord,
  K_fixed = 5,
  platform = "ST"
)

c_mat <- matrix(0, nrow = n_iter, ncol = numOfdata)
for (i in 1:n_iter) {
  c_mat[i, ] <- result$clIds_mcmc[[i]]
}
C_new <- apply(c_mat, 2, function(x) {
  as.numeric(names(table(x)) [which.max(table(x))])
})

C_new

```


### 2. Human lymph node A1 data

The following code shows an example (the first real dataset study in the manuscript) that runs the function "BSMIND" in our package.

Import the required R package.


```R
library(BSMIND)
```

Load the example data stored in this package. The example data includes:

* data_1: Gene expression data processed by standardization and principal component analysis (PCA).
* data_2: Pre-filtered proteomic data standardized via centered log-ratio (CLR) transformation, which conforms to a normal distribution.
* coord: Spatial coordinates corresponding to manually annotated tissue domain identification.
* true_label: The manual annotation of the domain identification.


```R
A1_data <- load_A1_data()
data_1 <- A1_data$gene_pca
data_2 <- A1_data$ADT_normal
coord <- A1_data$coord
true_label <- A1_data$true_label
```

Run "BSMIND" function to obtained consensus labels.
Obtain consensus clustering labels (C_new).
```R
set.seed(73597)
result <- BSMIND(
  data_1 = data_1,
  data_2 = data_2,
  coord = coord,
  K_fixed = 10,
  platform = "Visium"
)

c_mat <- matrix(0, nrow = n_iter, ncol = numOfdata)
for (i in 1:n_iter) {
  c_mat[i, ] <- result$clIds_mcmc[[i]]
}
C_new <- apply(c_mat, 2, function(x) {
  as.numeric(names(table(x)) [which.max(table(x))])
})

C_new

```


 
## Contact

If you have any questions regarding this package, please contact Mingxiu Wei at [yinqiaoyan@bjut.edu.cn](mailto:weimingxiu@emails.bjut.edu.cn).































































