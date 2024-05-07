

This vignette provides all reproducible codes for our article:

# LUPINE

Saritha Kodikara and Kim-Anh Lê Cao 

**Abstract:** The microbiome is a complex ecosystem of interdependent taxa that has traditionally been studied through cross-sectional studies. Longitudinal microbiome studies are becoming increasingly popular. These studies enable researchers to infer taxa associations towards the understanding of coexistence, competition, and collaboration between bacteria across time. 
However, traditional metrics for association analysis, such as correlation, are limited due to the data characteristics of microbiome data (sparse, compositional, multivariate). Several network inference methods have been proposed, but have been largely unexplored in a longitudinal setting.

We introduce LUPINE (LongitUdinal modelling with Partial least squares regression for NEtwork inference), a novel approach that leverages on conditional independence and low-dimensional data representation. This method is specifically designed to handle scenarios with small sample sizes and small number of time points. LUPINE generates networks while considering information from all past time points. It is the first method of its kind to infer microbial networks across time. We validate LUPINE and its variant, LUPINE\_single (for a single time point) in simulated data and three case studies, where we highlight LUPINE's ability to identify relevant taxa in each study context, across different experimental designs (mouse and human studies, with and without interventions, short and long time courses). To detect changes in the networks across time or in response to external disturbances, we used different metrics to compare our inferred networks.

LUPINE is a simple yet innovative network inference methodology that is suitable for, but not limited to, analysing longitudinal microbiome data. The R code and data are publicly available for readers interested in applying these new methods in their studies.



**Keywords:**  longitudinal, microbiome networks, partial correlation, 16S
