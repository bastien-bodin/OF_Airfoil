# Automated NACA Airfoil CFD: Python & OpenFOAM

**Author:** Bastien Bodin (PhD)  
**GitHub:** [bastien-bodin](https://github.com/bastien-bodin)  
**Contact:** bastien.bodin@proton.me

---

## 📖 Overview

This project provides a **fully automated, open-source workflow** to perform Computational Fluid Dynamics (CFD) simulations on **NACA 4-digit airfoils**.

Unlike traditional workflows relying on commercial meshing tools (Pointwise, Ansys ICEM), this project utilizes **Python** and the **Gmsh API** to generate high-quality, structured-like **C-Grid meshes**. The pipeline automatically handles the conversion to OpenFOAM, patches the boundary conditions for 2D simulations, and sets up the numerical schemes for a steady-state RANS analysis (`simpleFoam`).

---

## 🚀 Key Features

### Analytical Geometry Generation
- Exact mathematical reconstruction of **NACA 4-digit profiles** (e.g., NACA 2412, 0012).

### Cosine Spacing
- Cluster points at the **leading and trailing edges** to capture high gradients.

### Automated Meshing (Gmsh)
- Generates a **C-Grid topology** ideal for external aerodynamics.
- Handles the **2D-to-3D extrusion** required by OpenFOAM (1-cell thick).

### OpenFOAM Bridge
- **Auto-correction** of boundary files (patches `frontAndBack` to type `empty` for pure 2D).
- Automatic generation of `controlDict`, `fvSchemes`, and `fvSolution` optimized for incompressible RANS.

### Commercial-Free
- Relies **100%** on open-source tools.

---

## 📐 Physics & Domain Strategy

### Domain Size Justification

The fluid domain dimensions are chosen based on the **asymptotic behavior of the velocity potential** in subsonic flows. Unlike thermal diffusion (which decays exponentially), incompressible pressure fields are **elliptical** and decay **algebraically**.

- **Upstream/Farfield** ($R \approx 20c$): To minimize blockage effects. The lift-induced velocity perturbation decays as $1/r$ (vortex singularity). A boundary placed too close enforces $U = U_\infty$ artificially, violating the circulation requirement.

- **Downstream** ($L \approx 30c$): Extended to allow the turbulent wake to diffuse before reaching the outlet patch, preventing reverse flow and numerical instability.

### Numerical Schemes

- **Time:** Steady-state (`simpleFoam`).
- **Convection:** Second-order `linearUpwind` for velocity (crucial for accurate lift/drag prediction), first-order `upwind` for turbulence quantities (stability).
- **Pressure-Velocity Coupling:** SIMPLEC algorithm (`consistent yes`) for faster convergence.

---

## 🛠️ Prerequisites

Ensure you have the following installed on your **Linux/WSL** system:

1. **OpenFOAM** (Tested on v2412, compatible with Foundation/ESI versions).
2. **Python 3.x** with the following libraries:
   ```bash
   pip install numpy gmsh matplotlib
   ```
3. **System Libraries for Gmsh:**  
   If you encounter `libGLU` errors:
   ```bash
   sudo apt-get install libglu1-mesa  # Ubuntu/Debian
   ```

---

## ⚡ Usage

### Clone the repository
```bash
git clone https://github.com/bastien-bodin/openfoam-airfoil-python.git
cd openfoam-airfoil-python
```

### Run the generation script
This script generates the mesh, converts it, and writes the OpenFOAM configuration files.
```bash
python3 main.py
```
*By default, it generates a NACA 2412 mesh.*

### Run the Simulation
Load your OpenFOAM environment and launch the solver:
```bash
# Source OpenFOAM (if not in .bashrc)
source /opt/openfoam2412/etc/bashrc

# Run the solver
simpleFoam
```

### Post-Processing
Visualize the results using ParaView:
```bash
paraFoam
```

---

## 📂 Directory Structure

```
.
├── main.py                 # Master script (Geometry + Mesh + Case Setup)
├── airfoil.msh             # Generated Gmsh mesh (artifact)
├── system/
│   ├── controlDict         # Time and I/O control (Auto-generated)
│   ├── fvSchemes           # Discretization schemes (Auto-generated)
│   └── fvSolution          # Linear solver settings (Auto-generated)
├── constant/
│   └── polyMesh/           # Mesh files (boundary, faces, points...)
└── 0/                      # Initial conditions (U, p, nut, k, omega...)
```

---

## 📝 To-Do / Roadmap

- [ ] Implement k-omega SST turbulence model setup automatically.
- [ ] Add automatic calculation of Lift ($C_L$) and Drag ($C_D$) coefficients.
- [ ] Parameterize Angle of Attack (AoA) directly in the Python script.

---

## 📜 License

This project is licensed under the **MIT License**.