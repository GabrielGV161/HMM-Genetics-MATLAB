# HMM-Genetics-MATLAB
Hidden Markov Models for SNP Genetic Sequence Analysis
# Hidden Markov Models for Genetic Sequence Analysis

## Overview
This project explores the use of Hidden Markov Models (HMMs) to probabilistically identify genomic regions (SNPs) associated with a phenotypic trait. The work focuses on modeling dependencies between consecutive SNPs, leveraging Markovian structure to improve genetic mapping beyond independent SNP analysis.

The project was developed as an academic MATLAB assignment and later extended and structured for reproducibility and clarity.

## Objectives
- Model SNP–trait associations using a Hidden Markov framework
- Implement and analyze:
  - Forward algorithm (probability accumulation)
  - Viterbi algorithm (most likely hidden state sequence)
- Study how transition probabilities affect genomic region inference
- Visualize probabilistic state evolution along SNP sequences

## Methods
- Definition of hidden states and emission probabilities
- Implementation of forward-backward and Viterbi algorithms
- Parameter estimation using MATLAB built-in functions
- Validation using simulated and example genetic sequences
## Model Description
Hidden states (3):
- SNP associated with the progenitor carrying the trait
- SNP associated with the progenitor without the trait
- SNP not associated (random inheritance)
Observations:
- Simulated SNP indices representing genomic positions
Transition & emission probabilities:
- Defined heuristically to reflect biological intuition (linkage disequilibrium).

Note: No Expectation-Maximization (EM) was used due to the absence of real genomic datasets.


This work prioritizes conceptual understanding and modeling design rather than parameter estimation from large biological cohorts.

## Tools
- MATLAB
- Probabilistic modeling
- Genetics / bioinformatics concepts

## Implementation details
Custom MATLAB scripts implementing:
- A custom function, prrMINE2, was designed and implemented to adapt the classical Forward Algorithm to the specific structure and assumptions of the proposed genomic HMM.
- Viterbi decoding
- MATLAB built-in HMM functions (hmmgenerate, hmmviterbi) used only for sequence generation and validation, not as a black- box solution
## Custom visualization:
- Heatmaps
- Probability evolution plots
- State comparison under different transition matrices

## Results
The model successfully computes state probability distributions across genomic sequences and provides insight into the most likely hidden state paths.

Results demonstrate how probabilistic modeling can be applied to biological sequence data even under limited data availability, while highlighting the importance of parameter estimation methods for future improvements.

### Current Limitations
- Transition probabilities were defined heuristically
- No parameter learning was performed due to lack of sufficient data

## Author
Gabriel Garrido Vallet
## Authorship:
This project was developed as part of a group academic assignment.
This repository reflects my personal contribution to the group academic project.
The original work was developed collaboratively, with my main contributions
covering model design, implementation, analysis, and documentation.

## Date
2024
