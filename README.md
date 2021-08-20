# MuGLasso_GWAS
This repository is an implementation of MuGLasso applied to GWAS, the method is presented in the paper [Multitask group Lasso for Genome Wide association Studies in admixed populations](https://www.biorxiv.org/content/10.1101/2021.08.02.454499v1.full.pdf).
## Simulated data 
We used [GWAsimulator](https://biostat.app.vumc.org/GWAsimulator) tool to generate simulated GWAS data of admixed populations (CEU and YRI) in plink format.
Here are the following steps: 
- Prepare the input control files, we used ```Files/control_CEU.txt``` and ```Files/control_YRI.txt``` to obtain the used dataset in the paper.
- Prepare the input reference data in phased format (to specify in the first line in the controls file). In our case, we used HapMap3 for both subpopulations (CEU and YRI) as input reference for GWAsimulator, it could be downloaded in: https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII/.
You can use other reference dataset of your choice.
- Generate two subsets for each subpopulation: 
```
GWAsimulator control_ceu.txt [seed number]
GWAsimulator control_yri.txt [seed number]
```
The obtained output data is in plink format, we merge both populations using [plink](https://www.cog-genomics.org/plink/) software (with ```--bmerge``` flag). Please make sure that you change the Individual IDs for one subpopulation, as both subpopulations files contains the same IDs, this can be handled in plink using the flag ```--update-ids```.
#### Quality control and preprocessing 

Quality control steps are required to remove poorly performing SNPs, the following filters were used to obtain the data in the paper: 
```
plink --bfile [data] --maf 0.05 --hwe 1e-4 --make-bed --out [data_qc]
```

We also did Linkage Disequilibrium pruning to remove strong correlation between SNPs: 
```
plink --bfile [data_qc] --indep-pairwise 50 5 0.85 --out tmp
plink --bfile [data_qc] --extract tmp.prune.in --make-bed --out [data_qc_ld]
```

More details about how to use GWAsimulator in the github repository: [GWAS-admixed-population-simulator](https://github.com/asmanouira/GWAS-admixed-population-simulator).

## Linkage Disequilibrium groups clustering 

Hierarchical clustering was used to partitionnate the SNPs in blocks following LD patterns. 
To obtain the LD-groups for each subpopulation, we used [Adjclust](https://github.com/pneuvial/adjclust) software. 
The R code is given in ```adjclust.R``` 
We then combine those LD-groups into common shared LD-groups for both subpopulations using ```get_shared_groups.py```.

## MuGLasso_GWAS

**MuGLasso** is a multitask approach that performs feature selection at the LD-groups level using group Lasso. The input tasks correspond to subpopulations (*two tasks in our simulations case: for CEU and YRI*).

#### Architecture 
![Image description](/Images/MuGLasso.jpg)

#### Problem formulation and implementation
The problem is reformulated as a group Lasso with G groups each containing T copies (one per task) of the features of SNPs group.

The reformulated data is obtained in ```input_matrix.py```.

The model is implemented in ```MuGLasso.py``` and it is an adaptation of the sparse group Lasso of [Gap_safe_rules](https://github.com/EugeneNdiaye/Gap_Safe_Rules) package, it is a sparse group Lasso with squared loss. 

In our implementation, we extend it for a **group** Lasso with **logistic** loss as the phenotype in our dataset is categorical (case-control study). 
**Gap safe screening rules** make the solvers faster and avoid memory errors by discarding useless coefficients in the optimization procedure. More details in their paper: [Gap Safe screening rules for sparsity enforcing penalties](https://arxiv.org/abs/1611.05780).

To improve the stability of the model, we integrate subsampling following [Meinshausen procedure](https://rss.onlinelibrary.wiley.com/doi/full/10.1111/j.1467-9868.2010.00740.x). Code in ```stab_sel.py```.
#### Evaluation

**Pearson's stability index:** we estimate the stability of the selection of MuGLasso by computing [Pearson's stability](http://www.cs.man.ac.uk/~nogueirs/files/ecml2016.pdf) index implemented in ```get_stability.py  ```.



*Note that a tutorial will be made shortly available in a notebook.*
