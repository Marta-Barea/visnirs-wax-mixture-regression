# **Machine Learning for VIS-NIR Wax Mixture Regression**

The **visnirs-wax-mixture-regression** repository includes data processing methods and machine learning approaches applied to the analysis and characterization of petroleum wax blends. The focus is on classifying and predicting the composition of wax mixtures using VIS-NIR spectroscopy data.

---

## 🛠️ **System Requirements**

### Software
- **R version 4.1.2** (or higher)
- RStudio (optional but recommended)

1. **Feature Selection**:
   - Boruta Algorithm using the `Boruta` package (v7.0.0).
   - Genetic Algorithm implemented via the `caret` package (v6.0-90).

2. **Supervised Learning Models**:
   - Developed with the `caret` package (v6.0-90).
   - Includes Partial Least Squares Regression (PLSR), Support Vector Regression (SVR), and Random Forest Regression (RF).
   - Evaluation metrics calculated using the `MLmetrics` package (v1.1.1).

3. **Data Preprocessing**:
   - Spectral smoothing and normalization using the `prospectr` package (v0.2.3).

4. **Visualization**:
   - Results visualized with the `ggplot2` package (v3.3.5) and base R graphics.

5. **Interactive Web Application**:
   - Built with the `shiny` package (v1.7.1).

---

## 📂 Repository Structure

- `spectra/`: Contains scripts for raw and preprocessed VIS-NIR spectral data visualization.
- `feature selection plot/`: Scripts for visualizing feature importance from Boruta and Genetic Algorithm.
- `supervised algorithms/`: Scripts for supervised learning experiments.
- `unsupervised algorithms/`: Scripts for unsupervised learning and clustering analysis.
- `App/`: Shiny application for interactive exploration of the results.
- `LICENSE`: Licensing information.
- `README.md`: This file.

---
## ⚙️ How to Use This Repository

### Clone the Repository


   <pre markdown="1"> ```bash
git clone https://github.com/Marta-Barea/visnirs-wax-mixture-regression
cd visnirs-wax-mixture-regression </pre>

### Running the Shiny Application
1. Place `app.R`, `svm.rds`, `svr.rds`and `test_data.xlsx` in the same folder.
2. In your R console, run: 
   
   <pre markdown="1"> ```R 
      shiny::runApp("app.R") </pre>

3. Use the web interface to:
- 📁 **Upload** `.csv` or `.xlsx` data files.
- 🛠️ **Preprocess** data using advanced filtering techniques.
- 🤖 **Predict** hydroprocessing grades with AI.

---

### 📂 **Example Dataset**
A sample dataset (`test_data.xlsx`) is included for demonstration purposes. It contains Vis-NIR spectral readings and hydroprocessing grades for various wax samples.


---

### 🤝 **Contributors**
- **University of Cádiz (AGR-291 Research Group)**
  - Specializing in hydrocarbon characterization and spectroscopy.

---

### 📜 **License**
This project is licensed under the GNU GENERAL PUBLIC License. See `LICENSE` for details.
