# Machine learning-based approaches to Vis-NIR data for the automated characterization of petroleum wax blends

## Description

This repository contains the source code for all data processing and the application of machine learning algorithms used in the article "Machine learning-based approaches to Vis-NIR data for the automated characterization of petroleum wax blends".

## Contents

- `spectra/`: Folder containing the spectra data.
- `ffeature selection plot: Source code for visualizing variables selected by the Boruta Algorithm and Genetic Algorithm.
- `supervised algorithms/`: Source code for all the supervised machine learning models and experiments.
- `unsupervised algorithms/`: Source code related to unsupervised learning techniques and clustering.
- `App/`: A Shiny application to demonstrate and visualize the findings.

## Requirements

All data analysis was performed with **R (version 4.1.2)**. The software and packages used include:


- **prospectr (version 0.2.3)**: Implemented for spectral data processing using the SG algorithm.
- **Boruta (version 7.0.0)**: Employed for feature selection through the Boruta algorithm.
- **caret (version 6.0–90)**: Used for the application of the Genetic Algorithm (GA) and for the development of PLS, SVR, and RF models.
- **MLmetrics (version 1.1.1)**: Provided model evaluation metrics.
- **graphics (version 4.1.2)** and **ggplot2 (version 3.3.5)**: Facilitated data and model output visualization.
- **shiny (version 1.7.1)**: Enabled the development of an interactive web application.
