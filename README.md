# Nonlinear Cantilever Beam Analysis (Ritz Method + Newton-Raphson)

This repository contains a MATLAB implementation of a **nonlinear finite element Ritz solver** for a cantilever beam subjected to a tip load.  
The formulation uses **geometrically nonlinear kinematics (von Kármán strains)** and solves the nonlinear equilibrium equations using the **Newton–Raphson method**.

---

## 🔹 Features
- Ritz method with:
  - 5th-order polynomial approximation for **axial displacement (u)**
  - 3rd-order polynomial approximation for **transverse displacement (w)**
- Single-element, global Ritz approximation
- Explicit Newton–Raphson solver with load stepping
- Outputs:
  - **Load–deflection curve** (linear vs nonlinear)
  - **Final deformed beam shape** (linear vs nonlinear)

---

## 🔹 Theory
The total potential energy is formulated as:

\[
\Pi = \int_0^L \left[ \tfrac{1}{2} EA \, \varepsilon_0^2 + \tfrac{1}{2} EI \, \varepsilon_1^2 \right] dx - P w(L)
\]

where  
- \( \varepsilon_0 = u'(x) + \tfrac{1}{2}(w'(x))^2 \)  (axial strain)  
- \( \varepsilon_1 = -w''(x) \)  (curvature)  

The Ritz method substitutes polynomial trial functions into the energy functional.  
Equilibrium is obtained by enforcing stationarity via Newton–Raphson iterations.

---

## 🔹 Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/nonlinear-beam-ritz.git
   cd nonlinear-beam-ritz
