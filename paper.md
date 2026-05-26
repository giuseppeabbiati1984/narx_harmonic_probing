---
title: 'Numerical harmonic probing of nonlinear systems'
tags:
  - MATLAB
  - System identification
  - Nonlinear dynamics
  - Harmonic probing
  - NARX
authors:
  - name: David Stamenov
    orcid: 0009-0004-5718-9884
    corresponding: true  
    affiliation: 1
  - name: Thomas Sauder
    orcid: 0000-0001-7445-7239
    equal-contrib: true 
    affiliation: "2, 3"
  - name: Giuseppe Abbiati
    orcid: 0000-0002-5048-8505
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Aarhus University, Inge Lehmanns Gade 10, 8000, Aarhus, Denmark
   index: 1
   ror: https://ror.org/01aj84f44
 - name: SINTEF Ocean, P.O. Box 4762 Torgarden, 7465, Trondheim, Norway
   index: 2
   ror: https://ror.org/004wre089
 - name: Department of Marine Technology, Norwegian University of Science and Technology, Trondheim, Norway
   index: 3
   ror: https://ror.org/05xg72x27
date: 22 May 2026
bibliography: paper.bib

---

# Summary

Nonlinear engineering systems exhibit frequency intermodulation, whereby energy present in the input spectrum is transferred to additional frequency components in the system response. Such behavior cannot be represented by linear frequency response functions and instead requires higher-order frequency-domain descriptions, such as generalized frequency response functions derived from Volterra-series theory. The numerical harmonic probing algorithm provides a practical route for extracting these functions from nonlinear input–output models. This software presents a MATLAB code implementation for estimating generalized frequency response functions from polynomial NARX models trained on measured input–output data. The method is formulated as a recursive numerical procedure in which harmonic balance residuals are evaluated directly at selected probing frequencies.

# Statement of need

Engineering systems are often approximated by linear models, which can provide an adequate representation of the underlying dynamics over a limited operating range. A defining property of linear systems is that each frequency component present in the input appears in the output at the same frequency, modified only in amplitude and phase according to the linear frequency response function. Many natural and engineered systems, however, exhibit nonlinear behavior. In such systems, energy may be transferred between frequencies, producing response components at frequencies that are not present in the input. This phenomenon, commonly referred to as frequency intermodulation, is observed in a broad range of engineering applications, including hydrodynamic low-frequency and high-frequency loads \cite{sauder_empirical_2021, stansberg_cross-bi-spectral_1994}, geometrically nonlinear vibrations \cite{Lacarbonara_nonlinear_2007}, cracked beam-like structural elements \cite{Surace_detecting_2011}, nonlinear aeroelastic bridge response under turbulent wind excitation \cite{Chen_aeroelastic_2003, Kravrakov_comparative_2017}, among many other. \par

Since linear models do not predict intermodulated frequency components, they are insufficient for analyzing the above-mentioned natural phenomena. As a result, engineers and scientists have been focused on developing and improving methods for analyzing nonlinear systems. \par

An important nonlinear system identification approach relies on the functional series of Volterra and Wiener \cite{billings_identification_1980, schetzen_1980_volterra}. The theoretical basis for these methods was founded by the pioneering works of Volterra, who introduced the theory of analytical functions \cite{volterra_sopra_1887} and later Frechet, who used it to represent continuous nonlinear systems \cite{frechet_sur_1910}. The early efforts of identifying nonlinear models, notably the works of Wiener, Bose, and Barrett \cite{wiener_nonlinear_1958, bose_theory_1956, barrett_hermite_1964}, revolved around a direct estimation of the Volterra kernels, also referred to as generalized frequency response functions (GFRFs), from measured data. However, estimating GFRFs directly from raw data requires a large number of identified parameters, namely kernel points, to adequately represent the Volterra kernels, and consequently demands large amounts of input–output measurement data. \par
%
A practical remedy to this is to identify the Volterra kernels indirectly by applying \emph{harmonic probing} \cite{bedrosian_output_1971} to data-driven mathematical input–output models. In this way, the kernels are not identified pointwise, but are instead characterized through a comparatively small set of trained model coefficients. Among the most widely used model classes for this purpose is the NARMAX model (nonlinear autoregressive moving average model with exogenous inputs) \cite{Billings_1985_partI, Billings_1985_partII, billings_nonlinear_2013}, together with the frequently adopted simpler NARX model \cite{chen_non-linear_1990}, in which the stochastic component is represented as additive noise acting only at the system output. \par
%
The harmonic probing method is traditionally carried out either through manual derivations or with the aid of symbolic toolboxes, both of which are subject to practical limitations in terms of scalability, efficiency, and ease of implementation. To address these limitations, the authors developed a numerical formulation of the harmonic probing method \cite{stamenov_numerical_2025} which provides the core theoretical foundation behind the software presented here. \par

The software enables efficient analysis of complex nonlinear systems for which developing an explicit physical model is impractical, but representative input–output data can be measured. The user provides measured input and output data, defines the NARX model settings, the operational frequency grid, and the software performs the complete workflow from data preparation to model fitting, numerical probing, and visualization of the resulting GFRFs. A preliminary implementation of the software has been applied and experimentally validated for the identification of nonlinear hydrodynamic loading models \cite{stamenov_comparison_2025}.

# State of the field                                                                                                                  

Several tools exist for galactic dynamics computations:                                                     
`galpy` [@Bovy:2015] is a Python package with similar goals,
providing orbit integration and potential classes for galactic dynamics.                                                              
`NEMO` [@Teuben:1995] is a well-established, comprehensive stellar dynamics                                                           
toolbox written primarily in C, offering extensive functionality but with a                                                           
steeper learning curve and less integration with modern Python workflows.                                                             
Other tools like `GalPot` provide specific Milky Way potential models but lack                                                        
the broader dynamical analysis capabilities.                                                                                          
                                                                                                                                        
`Gala` was built rather than contributing to existing projects for several                                                            
reasons. First, `Gala` was designed from the ground up to integrate seamlessly                                                        
with the Astropy ecosystem, using `astropy.units` and `astropy.coordinates`                                                           
as core dependencies rather than optional features. This tight integration                                                            
enables natural workflows for astronomers already using Astropy. Second,                                                              
`Gala`'s object-oriented API with consistent interfaces across subpackages                                                            
(potentials, integrators, dynamics) provides a more modular and extensible                                                            
design than alternatives available at the time. Third, `Gala` fills a specific                                                        
niche between simple demonstration codes and full N-body simulation packages                                                          
like `Gadget` [@Springel:2005] – it focuses on the common tasks in galactic                                                             
dynamics research (orbit integration, potential evaluation, coordinate                                                                
transformations) while maintaining both performance through C implementations                                                         
and usability through its Python interface.  

# Software design

\subsection{Software architecture}
The software is a MATLAB implementation for estimating generalized frequency response functions from time-series data using the numerical harmonic probing framework described in \cite{stamenov_numerical_2025}. The software is organized as sequential blocks starting from measured input-output data and ending with estimated GFRFs. The process consists of four main layers: data preparation, model fitting, numerical probing, and visualization. Figure~\ref{fig:method_sequence} shows the sequence of operations starting from a nonlinear system of interest and concluding with the GFRFs.
\begin{figure} [htb]
    \centering
    \includegraphics[width=1.1\linewidth]{Figures/method_sequence.png}
    \caption{Workflow from measured input–output data of a nonlinear system to its frequency-domain representation. The process includes data preparation, model fitting, numerical harmonic probing, and visualization of the resulting GFRFs.}
    \label{fig:method_sequence}
\end{figure} 
%
\subsection*{Data preparation}
In this stage, the measured input and output signals are prepared to ensure that they are suitable for nonlinear model fitting and subsequent probing. The raw time series data are resampled to a user-specified sampling rate and arranged into data segments of a user-prescribed length. To improve numerical conditioning and make the fitted coefficients comparable across datasets, the signals are normalized with respect to their standard deviations. The corresponding scaling factors are stored and later used to recover the transfer functions in physical units after probing.
%
\subsection*{Model fitting}
In the second stage, the measured input-output relationship is represented by a polynomial NARX model. The settings of the NARX model are defined by the number of lags on the input and output and the order of the model, all of which are input by the user. The polynomial NARX model is fitted by first assembling a regression matrix containing all candidate monomials constructed from the selected lagged input and output samples up to the prescribed polynomial order. Since this candidate set may become large and highly redundant, the coefficient vector is estimated using LASSO regularization \cite{tibshirani_regression_1996}. Specifically, the fitting problem is formulated as the minimization of the one-step-ahead prediction error augmented by an \(\ell_1\)-penalty on the model coefficients. This penalty promotes sparsity by shrinking weak or non-informative terms to zero, thereby performing coefficient estimation and regressor selection simultaneously. The resulting sparse polynomial NARX model retains only the dominant terms needed to describe the nonlinear input-output relationship, which improves robustness and interpretability and provides a compact basis for the subsequent numerical harmonic probing procedure.
\begin{figure}
    \centering
    \includegraphics[width=1\linewidth]{Figures/NARX_fit.png}
    \caption{Workflow for the fitting of the NARX model}
    \label{fig:NARX fit}
\end{figure}

\subsection*{Numerical harmonic probing}
Once the NARX model is available, the rest of the workflow is independent of the particular physical system that generated the data. The fitted polynomial NARX model is transformed from a time-domain representation into a frequency-domain description by recursively estimating the generalized frequency response functions. A discrete probing frequency grid is first defined, after which the GFRFs are computed sequentially. At each selected frequency combination, the software constructs a trial harmonic response, inserts it into the polynomial model, and evaluates the resulting residual at the target output frequency. Since this residual depends linearly on the unknown transfer function value, the latter can be recovered numerically from a small number of residual evaluations, without deriving symbolic harmonic balance expressions. In this way, the probing stage systematically maps the fitted nonlinear input-output model into a set of frequency response functions that quantify the dynamics of the system of interest.
%
\subsection*{GFRF visualization}
Lastly, the estimated GFRFs are organized and displayed to enable direct interpretation of the identified nonlinear dynamics. After the probing has been completed for all data segments, the software computes summary statistics, such as the mean and standard deviation of the estimated transfer functions across segments, in order to assess the consistency and robustness of the results. The quantities are then plotted over the user-defined probing frequency domain, by showing the identified response together with uncertainty bands. For higher-order transfer functions, selected diagonals or slices of the multidimensional frequency domain may be visualized to highlight dominant interaction patterns in a compact and interpretable form.
%
\begin{figure} [H]
    \centering
    \includegraphics[width=1\linewidth]{Figures/flowchart_HP.png}
    \caption{Workflow for the numerical harmonic probing up to third-order.}
    \label{fig:flowchart_HP}
\end{figure}

# Research impact statement

The main impact of the software lies in making higher-order frequency-domain analysis of nonlinear systems more accessible to the system identification community. While nonlinear black-box models such as NARX are flexible and can often be fitted successfully to measured input-output data, their coefficients are usually difficult to interpret physically. By converting a fitted nonlinear input-output model into generalized frequency response functions (GFRFs), the software provides a practical route from time-domain identification to an interpretable frequency-domain description. GFRFs reveal how nonlinear systems redistribute energy across frequencies through intermodulation and thereby offering insight that is not directly visible from the polynomial coefficients of the time-domain model alone. \par

The numerical probing workflow enables the identified model to be inspected through GFRFs which makes it easier to diagnose dominant nonlinear interactions and assess which frequency combinations are responsible for the observed response. In this sense, the software does not merely improve model estimation; it improves model interpretability, which is a central challenge in nonlinear system identification.

Because the probing is performed numerically and recursively, the framework is also well suited for parallel computing and future extension to higher-order or multi-scale systems. This is particularly important for the system identification community, where there is a growing need for tools that can handle complex nonlinear dynamics while remaining computationally tractable and sufficiently transparent for scientific analysis. Overall, the software contributes a practical and scalable pathway for transforming measured input-output data into physically interpretable nonlinear frequency-domain models.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References