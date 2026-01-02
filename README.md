# Reissner-Mindlin Plate Theory Calculator

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](YOUR_STREAMLIT_APP_URL_HERE)

A Python-based interactive tool for calculating deformations, stresses, and resultant forces in thick plates based on **Reissner-Mindlin Plate Theory**. This application digitizes the theoretical derivation from the *Aircraft Structures* course (University of Liège), providing a transparent "Constitutive Law Calculator" for students and engineers.

## 🚀 Live Demo
**[Click here to try the application live](YOUR_STREAMLIT_APP_URL_HERE)** *No installation required.*

## 📖 Overview
This tool is designed to validate the **Constitutive Laws** (Material Laws) of plate mechanics. Unlike standard FEA software, this script explicitly logs every step of the calculation (with LaTeX formulas and page references), making it ideal for verifying hand calculations or understanding the underlying theory.

### Features
* **Kinematics:** Calculate strain ($\epsilon$), curvature ($\kappa$), and transverse shear ($\gamma$) from displacement gradients.
* **Stress Recovery:** Compute through-thickness stress distributions $\sigma^{\alpha\beta}(z)$ and elongation $\lambda_h$.
* **Resultants:** Integrate stresses to obtain Membrane Forces ($N$), Bending Moments ($M$), and Shear Forces ($Q$).
* **Equilibrium Check:** Verify the balance between internal force gradients and external loads.

## 🛠️ Installation (Local)
If you prefer to run the code locally:

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/reissner-mindlin-calculator.git](https://github.com/YOUR_USERNAME/reissner-mindlin-calculator.git)
    cd reissner-mindlin-calculator
    ```

2.  **Install dependencies:**
    ```bash
    pip install -r requirements.txt
    ```

3.  **Run the application:**
    ```bash
    streamlit run app.py
    ```

## 📚 Theory & References
This implementation is based strictly on the derivation provided in:
> **Aircraft Structures - Plates - Reissner-Mindlin Theory** > *Ludovic Noels, University of Liège (2013-2014)*

Key derivation pages referenced in the code:
* **Page 31-33:** Deformation Modes (Membrane, Bending, Shear)
* **Page 40-41:** Stress Recovery & Through-thickness elongation
* **Page 43-47:** Resultant Stiffness Tensors ($\mathcal{H}_n, \mathcal{H}_m, \mathcal{H}_q$)
* **Page 51:** Equilibrium Equations

## 📄 License
MIT License