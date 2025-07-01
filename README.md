# Master's Thesis: Functional Regression and Forecasting of River Flow

This repository contains the materials and code for my Master's thesis:

**Functional Regression and Forecasting of River Flow in Environmental Time Series in Ticino: A Spatial-Temporal Approach with OASI Data**

---

## 📂 Repository Structure

The repository is organized into folders corresponding to the chapters of the thesis:

- **Chapter 2** – Data exploration and preprocessing
- **Chapter 3** – Baseline modeling approaches
- **Chapter 4** – Functional regression methodology
- **Chapter 5** – Evaluation and forecasting results
- **App** – Shiny application for visualization and prediction

Each folder contains:
- R scripts used for data analysis and modeling
- Supporting files such as figures, outputs, and documentation

---

## 🛠️ Technologies Used

- **R** (data cleaning, modeling, and visualization)
- **Shiny** (interactive web application)
- **LaTeX** (thesis document preparation)

---

## 📖 Abstract 
The Osservatorio Ambientale della Svizzera Italiana is a tool that promotes continuous environmental monitoring, with modern and flexible data management, using sensors installed around Ticino. Having real-time, complete environmental data is of growing importance nowadays; hence, developing models that can handle incomplete data and accurately predict environmental variables is fundamental.
Since the data collected by the sensors often suffer from gaps, there is a need for robust systems that can adequately impute missing data, particularly in the area of water resources.
This thesis addresses the challenge of imputing and predicting daily river flow using different models, from simpler methods to advanced techniques such as a functional data analysis framework enriched with spatial information. Specifically, we develop a two-basis functional regression model that uses B-spline representations of precipitation curves as predictors for river flow curves. We start from a simple model that includes just precipitation, and gradually incorporate additional covariates such as time, elevation, and slope. The methodology is validated on a dataset from the OASI platform that span over multiple years, for both imputations and predictions.
Results show that the functional regression approach often outperforms standard imputation techniques, such as linear interpolation and GAM-based methods.
A Shiny application is also developed to make the model interactive and accessible.
