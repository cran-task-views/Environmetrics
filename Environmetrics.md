---
name: Environmetrics
topic: Analysis of Ecological and Environmental Data
maintainer: Gavin Simpson
email: ucfagls@gmail.com
version: 2023-04-05
source: https://github.com/cran-task-views/Environmetrics/
---

#### Introduction

This Task View contains information about using R to analyse ecological
and environmental data.

The base version of R ships with a wide range of functions for use
within the field of environmetrics. This functionality is complemented
by a plethora of packages available via CRAN, which provide specialist
methods such as ordination & cluster analysis techniques. A brief
overview of the available packages is provided in this Task View,
grouped by topic or type of analysis. As a testament to the popularity
of R for the analysis of environmental and ecological data, a [special
volume](http://www.jstatsoft.org/v22/) of the *Journal of Statistical
Software* was produced in 2007.

Those interested in environmetrics should consult the `r view("Spatial")` view.
Complementary information is also available in the `r view("Cluster")`,
and `r view("SpatioTemporal")` task views.

If you have any comments or suggestions for additions or improvements,
then please contact the maintainer or submit an issue or pull request
in the GitHub repository linked above.

A list of available packages and functions is presented below, grouped
by analysis type.

#### General packages

These packages are general, having wide applicability to the
environmetrics field.

-   Package `r pkg("EnvStats")` is the successor to the
    S-PLUS module *EnvironmentalStats*, both by Steven Millard. A [user
    guide in the form of a
    book](https://dx.doi.org/10.1007/978-1-4614-8456-1) has recently been
    released.

#### Modelling species responses and other data

Analysing species response curves or modelling other data often involves
the fitting of standard statistical models to ecological data and
includes simple (multiple) regression, Generalized Linear Models (GLM),
extended regression (e.g. Generalized Least Squares \[GLS\]),
Generalized Additive Models (GAM), and mixed effects models, amongst
others.

-   The base installation of R provides `lm()` and `glm()` for fitting
    linear and generalized linear models, respectively.
-   Generalized least squares and linear and non-linear mixed effects
    models extend the simple regression model to account for clustering,
    heterogeneity and correlations within the sample of observations.
    Package `r pkg("nlme")` provides functions for fitting
    these models. The package is supported by Pinheiro & Bates (2000)
    *Mixed-effects Models in S and S-PLUS*, Springer, New York. An
    updated approach to mixed effects models, which also fits
    Generalized Linear Mixed Models (GLMM) and Generalized non-Linear
    Mixed Models (GNLMM) is provided by the `r pkg("lme4")`
    package, though this is currently beta software and does not yet
    allow correlations within the error structure.
-   Recommended package `r pkg("mgcv", priority = "core")`
    fits GAMs and Generalized Additive Mixed Models (GAMM) with
    automatic smoothness selection via generalized cross-validation. The
    author of `r pkg("mgcv")` has also written a companion
    monograph, Wood (2017) *Generalized Additive Models; An Introduction
    with R* Second Edition. Chapman Hall/CRC, which has an accompanying package
    `r pkg("gamair")`.
-   Alternatively, package `r pkg("gam")` provides an
    implementation of the S-PLUS function `gam()` that includes LOESS
    smooths.
-   Proportional odds models for ordinal responses can be fitted using
    `polr()` in the `r pkg("MASS", priority = "core")`
    package, of Bill Venables and Brian Ripley.
-   A negative binomial family for GLMs to model over-dispersion in
    count data is available in `r pkg("MASS")`.
-   Models for overdispersed counts and proportions
    -   Package `r pkg("pscl")` also contains several
        functions for dealing with over-dispersed count data. Poisson or
        negative binomial distributions are provided for both
        zero-inflated and hurdle models.
    -   `r pkg("aod")` provides a suite of functions to
        analyse overdispersed counts or proportions, plus utility
        functions to calculate e.g. AIC, AICc, Akaike weights.
-   Detecting change points and structural changes in parametric models
    is well catered for in the `r pkg("segmented")` package
    and the `r pkg("strucchange")` package respectively.
    `r pkg("segmented")` is discussed in an
    [R News article](http://CRAN.R-project.org/doc/Rnews/Rnews_2008-1.pdf)

#### Tree-based models

Tree-based models are being increasingly used in ecology, particularly
for their ability to fit flexible models to complex data sets and the
simple, intuitive output of the tree structure. Ensemble methods such as
bagging, boosting and random forests are advocated for improving
predictions from tree-based models and to provide information on
uncertainty in regression models or classifiers.

##### Univariate trees

Tree-structured models for regression, classification and survival
analysis, following the ideas in the CART book, are implemented in

-   recommended package `r pkg("rpart")`
-   `r pkg("party")` provides an implementation of
    conditional inference trees which embed tree-structured regression
    models into a well-defined theory of conditional inference
    procedures

##### Multivariate trees

Multivariate trees are available in

-   package `r pkg("party")` can also handle multivariate
    responses.

##### Ensembles of trees

Ensemble techniques for trees:

-   The Random Forest method of Breiman and Cutler is implemented in
    `r pkg("randomForest")`, providing classification and
    regression based on a forest of trees using random inputs
-   Package `r pkg("ipred")` provides functions for improved
    predictive models for classification, regression and survival
    problems.

Graphical tools for the visualization of trees are available in package
`r pkg("maptree")`.

Packages `r pkg("mda")` and `r pkg("earth")`
implement Multivariate Adaptive Regression Splines (MARS), a technique
which provides a more flexible, tree-based approach to regression than
the piecewise constant functions used in regression trees.

#### Ordination

R and add-on packages provide a wide range of ordination methods, many
of which are specialized techniques particularly suited to the analysis
of species data. The two main packages are
`r pkg("ade4", priority = "core")` and
`r pkg("vegan", priority = "core")`.
`r pkg("ade4")` derives from the traditions of the French
school of "Analyse des Donnees" and is based on the use of the duality
diagram. `r pkg("vegan")` follows the approach of Mark Hill,
Cajo ter Braak and others, though the implementation owes more to that
presented in Legendre & Legendre (1988) *Numerical Ecology, 2<sup>nd</sup>
English Edition*, Elsevier. Where the two packages provide duplicate
functionality, the user should choose whichever framework that best
suits their background.

-   Principal Components (PCA) is available via the `prcomp()` function.
    `rda()` (in package `r pkg("vegan")`), `pca()` (in
    package `r pkg("labdsv", priority = "core")`) and
    `dudi.pca()` (in package `r pkg("ade4")`), provide more
    ecologically-orientated implementations.
-   Redundancy Analysis (RDA) is available via `rda()` in
    `r pkg("vegan")` and `pcaiv()` in
    `r pkg("ade4")`.
-   Canonical Correspondence Analysis (CCA) is implemented in `cca()` in
    both `r pkg("vegan")` and `r pkg("ade4")`.
-   Detrended Correspondence Analysis (DCA) is implemented in
    `decorana()` in `r pkg("vegan")`.
-   Principal coordinates analysis (PCO) is implemented in `dudi.pco()`
    in `r pkg("ade4")`, `pco()` in
    `r pkg("labdsv")`, `pco()` in
    `r pkg("ecodist")`, and `cmdscale()` in package
    `r pkg("MASS")`.
-   Non-Metric multi-Dimensional Scaling (NMDS) is provided by
    `isoMDS()` in package `r pkg("MASS")` and `nmds()` in
    `r pkg("ecodist")`. `nmds()`, a wrapper function for
    `isoMDS()`, is also provided by package
    `r pkg("labdsv")`. `r pkg("vegan")` provides
    helper function `metaMDS()` for `isoMDS()`, implementing random
    starts of the algorithm and standardized scaling of the NMDS
    results. The approach adopted by `r pkg("vegan")` with
    `metaMDS()` is the recommended approach for ecological data.
-   Coinertia analysis is available via `coinertia()` and `mcoa()`, both
    in `r pkg("ade4")`.
-   Co-correspondence analysis to relate two ecological species data
    matrices is available in `r pkg("cocorresp")`.
-   Canonical Correlation Analysis (CCoA - not to be confused with CCA,
    above) is available in `cancor()` in standard package stats.
-   Procrustes rotation is available in `procrustes()` in
    `r pkg("vegan")` and `procuste()` in
    `r pkg("ade4")`, with both `r pkg("vegan")`
    and `r pkg("ade4")` providing functions to test the
    significance of the association between ordination configurations
    (as assessed by Procrustes rotation) using permutation/randomization
    and Monte Carlo methods.
-   Constrained Analysis of Principal Coordinates (CAP), implemented in
    `capscale()` in `r pkg("vegan")`, fits constrained
    ordination models similar to RDA and CCA but with any
    dissimilarity coefficient.
-   Constrained Quadratic Ordination (CQO; formerly known as Canonical
    Gaussian Ordination (CGO)) is a maximum likelihood estimation
    alternative to CCA fit by Quadratic Reduced Rank Vector GLMs.
    Constrained Additive Ordination (CAO) is a flexible alternative to
    CQO which uses Quadratic Reduced Rank Vector GAMs. These methods and
    more are provided in Thomas Yee's `r pkg("VGAM")`
    package.
-   Fuzzy set ordination (FSO), an alternative to CCA/RDA and CAP, is
    available in package `r pkg("fso")`.
    `r pkg("fso")` complements a recent paper on fuzzy sets
    in the journal *Ecology* by Dave Roberts (2008, Statistical analysis
    of multidimensional fuzzy set ordinations. *Ecology* **89(5)**,
    1246-1260).

#### Dissimilarity coefficients

Much ecological analysis proceeds from a matrix of dissimilarities
between samples. A large amount of effort has been expended formulating
a wide range of dissimilarity coefficients suitable for ecological data.
A selection of the more useful coefficients are available in R and
various contributed packages.

Standard functions that produce, square, symmetric matrices of pair-wise
dissimilarities include:

-   `dist()` in standard package stats
-   `daisy()` in recommended package
    `r pkg("cluster", priority = "core")`
-   `vegdist()` in `r pkg("vegan")`
-   `dsvdis()` in `r pkg("labdsv")`
-   `Dist()` in `r pkg("amap")`
-   `distance()` in `r pkg("ecodist")`
-   a suite of functions in `r pkg("ade4")`

Function `distance()` in package `r pkg("analogue")` can be
used to calculate dissimilarity between samples of one matrix and those
of a second matrix. The same function can be used to produce pair-wise
dissimilarity matrices, though the other functions listed above are
faster. `distance()` can also be used to generate matrices based on
Gower's coefficient for mixed data (mixtures of binary, ordinal/nominal
and continuous variables). Function `daisy()` in package
`r pkg("cluster")` provides a faster implementation of
Gower's coefficient for mixed-mode data than `distance()` if a standard
dissimilarity matrix is required. Function `gowdis()` in package
`r pkg("FD")` also computes Gower's coefficient and
implements extensions to ordinal variables.

#### Cluster analysis

Cluster analysis aims to identify groups of samples within multivariate
data sets. A large range of approaches to this problem have been
suggested, but the main techniques are hierarchical cluster analysis,
partitioning methods, such as *k* -means, and finite mixture models or
model-based clustering. In the machine learning literature, cluster
analysis is an unsupervised learning problem.

The `r view("Cluster")` task view provides a more detailed
discussion of available cluster analysis methods and appropriate R
functions and packages.

Hierarchical cluster analysis:

-   `hclust()` in standard package stats
-   Recommended package `r pkg("cluster")` provides
    functions for cluster analysis following the methods described in
    Kaufman and Rousseeuw (1990) *Finding Groups in data: an
    introduction to cluster analysis*, Wiley, New York
-   `hcluster()` in `r pkg("amap")`
-   `r pkg("pvclust")` is a package for assessing the
    uncertainty in hierarchical cluster analysis. It provides
    approximately unbiased *p* -values as well as bootstrap *p* -values.

Partitioning methods:

-   `kmeans()` in stats provides *k* -means clustering
-   `cmeans()` in `r pkg("e1071")` implements a fuzzy
    version of the *k* -means algorithm
-   Recommended package `r pkg("cluster")` also provides
    functions for various partitioning methodologies.

Mixture models and model-based cluster analysis:

-   `r pkg("mclust")` and `r pkg("flexmix")`
    provide implementations of model-based cluster analysis.
-   `r pkg("prabclus")` clusters a species presence-absence
    matrix object by calculating an MDS from the distances, and applying
    maximum likelihood Gaussian mixtures clustering to the MDS points.
    The maintainer's, Christian Hennig, website contains several
    publications in ecological contexts that use
    `r pkg("prabclus")`, especially Hausdorf & Hennig (2007;
    [Oikos 116 (2007),
    818-828](https://doi.org/10.1111/j.0030-1299.2007.15661.x)).

#### Ecological theory

There is a growing number of packages and books that focus on the use of
R for theoretical ecological models.

-   `r pkg("vegan")` provides a wide range of functions
    related to ecological theory, such as diversity indices (including
    the "so-called" Hill's numbers \[e.g. Hill's N ^2^ \] and
    rarefaction), ranked abundance diagrams, Fisher's log series,
    Broken Stick model, Hubbell's abundance model, amongst others.
-   `r pkg("untb")` provides a collection of utilities for
    biodiversity data, including the simulation ecological drift under
    Hubbell's Unified Neutral Theory of Biodiversity, and the
    calculation of various diagnostics such as Preston curves.
-   Package `r pkg("BiodiversityR")` provides a GUI for
    biodiversity and community ecology analysis.
-   Function `betadiver()` in `r pkg("vegan")` implements
    all the diversity indices reviewed in Koleff et al (2003;
    [Journal of Animal Ecology 72(3),
    367-382](https://dx.doi.org/10.1046/j.1365-2656.2003.00710.x) ).
    `betadiver()` also provides a `plot` method to produce the
    co-occurrence frequency triangle plots of the type found in Koleff
    et al (2003).
-   Function `betadisper()`, also in `r pkg("vegan")`,
    implements Marti Anderson's distance-based test for homogeneity of
    multivariate dispersions (PERMDISP, PERMDISP2), a multivariate
    analogue of Levene's test (Anderson 2006; [Biometrics 62,
    245-253](https://dx.doi.org/10.1111/j.1541-0420.2005.00440.x) ).
    Anderson et al (2006; [Ecology Letters 9(6),
    683-693](https://dx.doi.org/10.1111/j.1461-0248.2006.00926.x) )
    demonstrate the use of this approach for measuring beta diversity.
-   The `r pkg("FD")` package computes several measures of
    functional diversity indices from multiple traits.

#### Population dynamics

##### Estimating animal abundance and related parameters

This section concerns estimation of population parameters (population
size, density, survival probability, site occupancy etc.) by methods
that allow for incomplete detection. Many of these methods use data on
marked animals, variously called 'capture-recapture',
'mark-recapture' or 'capture-mark-recapture' data.

-   `r pkg("Rcapture")` fits loglinear models to estimate
    population size and survival rate from capture-recapture data as
    described by [Baillargeon and
    Rivest (2007)](http://www.jstatsoft.org/v19/i05).
-   `r pkg("secr")` estimates population density given
    spatially explicit capture-recapture data from traps, passive DNA
    sampling, automatic cameras, sound recorders etc. Models are fitted
    by maximum likelihood. The detection function may be half-normal,
    exponential, cumulative gamma etc. Density surfaces may be fitted.
    Covariates of density and detection parameters are specified via
    formulae.
-   `r pkg("unmarked")` fits hierarchical models of
    occurrence and abundance to data collected on species subject to
    imperfect detection. Examples include single- and multi-season
    occupancy models, binomial mixture models, and hierarchical distance
    sampling models. The data can arise from survey methods such
    temporally replicated counts, removal sampling, double-observer
    sampling, and distance sampling. Parameters governing the state and
    observation processes can be modelled as functions of covariates.
-   Package `r pkg("RMark")` provides a formula-based R
    interface for the MARK package which fits a wide variety of
    capture-recapture models. See the [RMark
    website](http://www.phidot.org/software/mark/rmark/) and a [NOAA
    report](http://www.afsc.noaa.gov/Publications/ProcRpt/PR2013-01.pdf)
    (PDF) for further details.
-   Package `r pkg("marked")` provides a framework for
    handling data and analysis for mark-recapture.
    `r pkg("marked")` can fit Cormack-Jolly-Seber (CJS) and
    Jolly-Seber (JS) models via maximum likelihood and the CJS model
    via MCMC. Maximum likelihood estimates for the CJS model can be
    obtained using R or via a link to the Automatic Differentiation
    Model Builder software. A [description of the
    package](https://dx.doi.org/10.1111/2041-210X.12065) was published
    in Methods in Ecology and Evolution.
-   `r pkg("mrds")` fits detection functions to point and
    line transect distance sampling survey data (for both single and
    double observer surveys). Abundance can be estimated using
    Horvitz-Thompson-type estimators.
-   `r pkg("Distance")` is a simpler interface to
    `r pkg("mrds")` for single observer distance sampling
    surveys.
-   `r pkg("dsm")` fits density surface models to
    spatially-referenced distance sampling data. Count data are
    corrected using detection function models fitted using
    `r pkg("mrds")` or `r pkg("Distance")`.
    Spatial models are constructed as in `r pkg("mgcv")`.

Packages `r pkg("secr")` can also be used to simulate data
from the respective models.

See also the `r view("SpatioTemporal")` task view for analysis of animal
tracking data under *Moving objects, trajectories*.

##### Modelling population growth rates:

-   Package `r pkg("popbio")` can be used to construct and
    analyse age- or stage-specific matrix population models.

#### Environmental time series

-   Time series objects in R are created using the `ts()` function,
    though see `r pkg("tseries")` or
    `r pkg("zoo")` below for alternatives.
-   Classical time series functionality is provided by the `ar()`, and
    `arima()` functions in standard package stats for autoregressive
    (AR), moving average (MA), autoregressive moving average (ARMA) and
    integrated ARMA (ARIMA) models.
-   The `r pkg("forecast")` package provides methods and
    tools for displaying and analysing univariate time series forecasts
    including exponential smoothing via state space models and automatic
    ARIMA modelling
-   The `r pkg("dse")` package provide a variety of more
    advanced estimation methods and multivariate time series analysis.
-   Packages `r pkg("tseries")` and
    `r pkg("zoo")` provide general handling and analysis of
    time series data.
-   Irregular time series can be handled using package
    `r pkg("zoo")` as well as by `irts()` in package
    `r pkg("tseries")`.
-   `r pkg("pastecs")` provides functions specifically
    tailored for the analysis of space-time ecological series.
-   `r pkg("strucchange")` allows for testing, dating and
    monitoring of structural change in linear regression relationships.
-   Detecting change points in time series data &mdash; see
    `r pkg("segmented")` above.
-   The `r pkg("surveillance")` package implements
    statistical methods for the modelling of and change-point detection
    in time series of counts, proportions and categorical data. Focus is
    on outbreak detection in count data time series.
-   Package `r pkg("dynlm")` provides a convenient interface
    to fitting time series regressions via ordinary least squares
-   Package `r pkg("dyn")` provides a different approach to
    that of `r pkg("dynlm")`, which allows time series data
    to be used with any regression function written in the style of `lm`
    such as `lm()`, `glm()`, `loess()`, `rlm()` and `lqs()` from
    `r pkg("MASS")`, `randomForest()` (package
    `r pkg("randomForest")`), `rq()` (package
    `r pkg("quantreg")`) amongst others, whilst preserving
    the time series information.
-   The `r pkg("openair")` provides numerous tools to
    analyse, interpret and understand air pollution time series data
-   The `r pkg("bReeze")` package is a collection of widely
    used methods to analyse, visualize, and interpret wind data. Wind
    resource analyses can subsequently be combined with characteristics
    of wind turbines to estimate the potential energy production.
-   The `r pkg("Rbeast")` package provides a Bayesian model averaging method 
    to decompose time series into abrupt changes, trend, and seasonality and 
    can be used for changepoint detection, time series decomposition, and 
	nonlinear trend analysis.	

Additionally, a fuller description of available packages for time series
analysis can be found in the `r view("TimeSeries")` task
view.

#### Spatial data analysis

See the `r view("Spatial")` CRAN Task View for an overview
of spatial analysis in R.

#### [Extreme values]{#extremes}

`r pkg("ismev")` provides functions for models for extreme
value statistics and is support software for Coles (2001) *An
Introduction to Statistical Modelling of Extreme Values* , Springer, New
York. Other packages for extreme value theory include

-   `r pkg("evir")`
-   `r pkg("evd")`
-   `r pkg("evdbayes")`, which provides a Bayesian approach
    to extreme value theory
-   `r pkg("extRemes")`

See also the `r view("ExtremeValue")` task view for further information.


#### Phylogenetics and evolution

Packages specifically tailored for the analysis of phylogenetic and
evolutionary data include:

-   `r pkg("ape")`
-   `r pkg("ouch")`

UseRs may also be interested in Paradis (2006) *Analysis of
Phylogenetics and Evolution with R*, Springer, New York, a book in the
["Use R!" book series](https://www.springer.com/series/6991) from Springer.

#### Soil science

Several packages are now available that implement R functions for
widely-used methods and approaches in pedology.

-   `r pkg("soiltexture")` provides functions for soil
    texture plot, classification and transformation.
-   `r pkg("aqp")` contains a collection of algorithms
    related to modelling of soil resources, soil classification, soil
    profile aggregation, and visualization.
-   The Soil Water project on R-Forge.R-project.org provides packages
    providing soil water retention functions, soil hydraulic
    conductivity functions and pedotransfer functions to estimate their
    parameter from easily available soil properties. Two packages form
    the project: `r rforge("soilwaterfun")` and `r rforge("soilwaterptf")`.

#### Hydrology and Oceanography

A growing number of packages are available that implement methods
specifically related to the fields of hydrology and oceanography. Also
see the [Extreme Value](#extremes) and the [Climatology](#secclimatology)
sections for related packages.

-   `r pkg("hydroTSM")` is a package for management,
    analysis, interpolation and plotting of time series used in
    hydrology and related environmental sciences.
-   `r pkg("hydroGOF")` is a package implementing both
    statistical and graphical goodness-of-fit measures between observed
    and simulated values, mainly oriented to be used during the
    calibration, validation, and application of
    hydrological/environmental models. A related package is
    `r pkg("qualV")` which
    provides quantitative and qualitative criteria to compare models
    with data and to measure similarity of patterns
-   `r pkg("topmodel")` is a set of hydrological functions
    including an R implementation of the hydrological model TOPMODEL,
    which is based on the 1995 FORTRAN version by Keith Beven. New
    functionality is being developed as part of the
    `r rforge("RHydro")` package on R-Forge.
-   Package `r pkg("seacarb")` provides functions for
    calculating parameters of the seawater carbonate system.
-   Stephen Sefick's `r pkg("StreamMetabolism")` package
    contains function for calculating stream metabolism characteristics,
    such as GPP, NDM, and R, from single station diurnal Oxygen curves.
-   Package `r pkg("oce")` supports the analysis of
    Oceanographic data, including ADP measurements, CTD measurements,
    sectional data, sea-level time series, and coastline files.
-   The `r pkg("nsRFA")` package provides collection of
    statistical tools for objective (non-supervised) applications of the
    Regional Frequency Analysis methods in hydrology.
-   The `r pkg("boussinesq")` package is a collection of
    functions implementing the one-dimensional Boussinesq Equation
    (ground-water).
-   `r pkg("rtop")` is a package for geostatistical
    interpolation of data with irregular spatial support such as runoff
    related data or data from administrative units.

#### [Climatology]{#secclimatology}

Several packages related to the field of climatology.

-   `r pkg("seas")` implements a number of functions for
    analysis and graphics of seasonal data.
-   `r pkg("RMAWGEN")` is set of S3 and S4 functions for
    spatial multi-site stochastic generation of daily time series of
    temperature and precipitation making use of Vector Autoregressive
    Models.

#### Palaeoecology and stratigraphic data

Several packages now provide specialist functionality for the import,
analysis, and plotting of palaeoecological data.

-   Transfer function models including weighted averaging (WA), modern
    analogue technique (MAT), Locally-weighted WA, & maximum likelihood
    (aka Gaussian logistic) regression (GLR) are provided by the
    `r pkg("rioja")` and `r pkg("analogue")`
    packages.
-   Import of common, legacy, palaeo-data formats is provided by package
    `r pkg("vegan")` (Cornell format).
-   Stratigraphic data plots can be drawn using `Stratiplot()` function
    in `r pkg("analogue")` and functions `strat.plot()` and
    `strat.plot.simple` in the `r pkg("rioja")` package.
    Also see the `r github("paleolimbot/tidypaleo")`
    package, which provides tools to produce stratigraphic plots using
    `ggplot()`. A [blog
    post](https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
    by the maintainer of the
    `r github("paleolimbot/tidypaleo")` package, Dewey
    Dunnington, shows how to use the package to create stratigraphic
    plots.
-   `r pkg("analogue")` provides extensive support for
    developing and interpreting MAT transfer function models, including
    ROC curve analysis. Summary of stratigraphic data is supported via
    principal curves in the `prcurve()` function.

#### Other packages

Several other relevant contributed packages for R are available that do
not fit under nice headings.

-   Andrew Robinson's `r pkg("equivalence")` package
    provides some statistical tests and graphics for assessing tests of
    equivalence. Such tests have similarity as the alternative
    hypothesis instead of the null. The package contains functions to
    perform two one-sided t-tests (TOST) and paired t-tests of
    equivalence.
-   Thomas Petzoldt's `r pkg("simecol")` package provides
    an object oriented framework and tools to simulate ecological (and
    other) dynamic systems within R. See the [simecol
    website](https://www.simecol.de/) and a [R
    News](https://CRAN.R-project.org/doc/Rnews)
    [article](https://CRAN.R-project.org/doc/Rnews/Rnews_2003-3.pdf) on
    the package for further information.
-   Functions for circular statistics are found in
    `r pkg("CircStats")` and
    `r pkg("circular")`.
-   Package `r pkg("e1071")` provides functions for latent
    class analysis, short time Fourier transform, fuzzy clustering,
    support vector machines, shortest path computation, bagged
    clustering, naive Bayes classifier, and more\...
-   Package `r pkg("pgirmess")` provides a suite of
    miscellaneous functions for data analysis in ecology.
-   `r pkg("mefa")` provides functions for handling and
    reporting on multivariate count data in ecology and biogeography.
-   Sensitivity analysis of models is provided by package
    `r pkg("sensitivity")`.
    `r pkg("sensitivity")` contains a collection of
    functions for factor screening and global sensitivity analysis of
    model output.
-   Functions to analyse coherence, boundary clumping, and turnover
    following the pattern-based metacommunity analysis of [Leibold and
    Mikkelson (2002)](https://dx.doi.org/10.1034/j.1600-0706.2002.970210.x)
    are provided in the `r pkg("metacom")` package.
-   Growth curve estimation via noncrossing and nonparametric regression
    quantiles is implemented in package
    `r pkg("quantregGrowth")`. A supporting paper is [Muggeo
    et al. (2013)](https://dx.doi.org/10.1007/s10651-012-0232-1).
-   The `r pkg("siplab")` package provides an R platform for
    experimenting with spatially explicit individual-based vegetation
    models. A supporting paper is [García, O.
    (2014)](http://www.mcfns.com/index.php/Journal/article/view/6_36).
-   `r pkg("PMCMRplus")` provides parametric and
    non-parametric many-to-one and all-pairs multiple comparison
    procedures for continuous or at least interval based variables. The
    package provides implementations of a wide range of tests involving
    pairwise multiple comparisons.



### Links
-   [The vegan development site on GitHub](https://vegandevs.github.io/vegan/index.html)
-   [Thomas Yee's VGAM package for R](http://www.stat.auckland.ac.nz/~yee/VGAM/)
-   [The cocorresp development site on GitHub](https://gavinsimpson.github.io/cocorresp/)
-   [Brodgar](http://www.brodgar.com/)
-   [More information on the ade4 package can be found on the ADE4 website](http://pbil.univ-lyon1.fr/ADE-4/)
-   [Mike Palmer's Ordination website](http://ordination.okstate.edu/)
-   [Thomas Petzoldt's page about ecological modelling with R](https://wwwpub.zih.tu-dresden.de/%7Epetzoldt/)
-   [The FLR project web page for Fisheries Library in R.](https://www.flr-project.org/)
