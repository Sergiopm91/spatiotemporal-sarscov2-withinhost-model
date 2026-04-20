# Spatiotemporal SARS-CoV-2 Within-Host Model

A spatiotemporal within-host PDE model of SARS-CoV-2 infection integrating viral dynamics, immune responses, cytokine regulation, humoral immunity, and thermoregulatory behavior, solved numerically in MATLAB.

## Overview

This repository contains the MATLAB implementation of a nineteen-variable reaction--diffusion model for within-host SARS-CoV-2 dynamics. The framework integrates:

- viral replication and clearance,
- cytokine regulation,
- innate immune dynamics,
- adaptive cellular immunity,
- humoral response,
- and thermoregulatory behavior.

The model is designed to reproduce the main temporal and spatial features of acute infection, including early viral amplification, rapid innate activation, delayed adaptive control, bounded inflammatory response, localized infection foci, antibody depletion zones, delayed cytotoxic coverage, and transient thermal gradients.

## Main Features

- Nineteen coupled state variables
- Two-dimensional spatial domain
- Reaction--diffusion immune modeling
- Delayed adaptive activation through sigmoidal functions
- Saturating pyrogenic forcing for thermoregulatory response
- Automated figure generation for manuscript-ready outputs
- Cross-sectional spatial profiles
- Clinical-timeline comparison figure
- Temporal convergence analysis
- Shape-parameter sensitivity analysis
- Export of figures and summary tables

## Mathematical Structure

The model includes:

- Cytokines: IL-12, IFN-$\gamma$, TNF-$\alpha$, IL-6, IL-10
- Dendritic cells: naïve and activated
- Macrophages: naïve and activated
- T-helper cells: naïve, Th1, Th2
- Cytotoxic T cells: naïve and effector
- Viral and humoral variables: virus, infected cells, antibodies, B cells
- Thermoregulatory variable: local tissue temperature

The governing equations are solved over a two-dimensional tissue domain using an explicit finite-difference implementation in the present script, together with supplementary routines for shape-parameter analysis inspired by meshless RBF methodology.

## Repository Structure

```text
.
├── maincode.m
├── figures/
├── README.md
└── LICENSE
