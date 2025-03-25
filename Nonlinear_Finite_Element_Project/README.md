# Nonlinear Finite Element Implementation of Thermoelastic Material Model

## Overview
This project presents a nonlinear finite element implementation for a fully coupled thermoelastic material model. The model incorporates temperature-dependent material properties and considers the intrinsic coupling between mechanical deformation and heat transfer. The solution is achieved using a Newton-Raphson iterative method, and the implementation is verified against analytical solutions and Abaqus simulations.

## Features
- Fully coupled thermoelastic finite element analysis (FEA)
- Temperature-dependent material properties
- Implicit time discretization
- 2D quadrilateral elements with shape functions and Jacobian transformation
- Numerical integration using Gauss quadrature
- Newton-Raphson method for solving nonlinear coupled equations
- Validation with analytical solutions and Abaqus simulations

## Prerequisites
- **MATLAB**: The implementation is written in MATLAB.
- **Gmsh**: Used for mesh generation.
- **Abaqus (optional)**: For verification and comparison.

## Project Structure
```
├── mainCouplefem.m       # Main MATLAB script for running simulations
├── elementRoutine.m      # Element-level calculations
├── materialRoutine.m     # Material property definitions
├── loadMesh.m            # Mesh loading function
├── square100Element.inp  # Example input mesh file
├── square400Element.inp  # Example finer mesh input file
├── squareCenterHoleFineMesh.inp  # Mesh file for square plate with hole
├── README.md             # Project documentation
```

## Running the Code
1. Place all MATLAB `.m` files and `.inp` mesh files in a single directory.
2. Open MATLAB and navigate to the project directory.
3. Run `mainCouplefem.m` by typing:
   ```matlab
   mainCouplefem
   ```

## Configurations
Modify the following parameters in `mainCouplefem.m` to customize the analysis:
- **Load Case Selection:**
  ```matlab
  LoadCase = 5; % Choose from 5, 6, or 7 for different thermal loads
  ```
- **Material Properties:**
  ```matlab
  temperatureDependent = 'yes'; % Set 'yes' for temperature-dependent properties
  ```
- **Steady-State vs Transient Analysis:**
  ```matlab
  steadyState = 'yes'; % Set 'no' for transient analysis
  ```
- **Mesh Selection:**
  ```matlab
  readFile = fopen('square100Element.inp', 'r');
  ```

## Testing and Validation
### 1. Patch Test
Run a simple patch test to check for correctness:
```matlab
testCase = 3;
eigenValueTest = 'no';
numericalTesting = 'no';
temperatureDependent = 'no';
steadyState = 'yes';
```

### 2. Comparison with Analytical Solution
Set `testCase` to 1, 2, or 4 and compare with exact solutions:
```matlab
testCase = 1;
```

### 3. Numerical Testing
Run internal consistency tests:
```matlab
testCase = 1;
numericalTesting = 'yes';
```

### 4. Eigenvalue Test
Verify the stability of the thermal element:
```matlab
eigenValueTest = 'yes';
```

## Results
Results are visualized using MATLAB’s contour plots, showing:
- Displacement magnitude (`u`)
- Temperature distribution (`T`)
- Stress fields

## Future Work
- Extend the model to thermoelastic-plastic material behavior
- Implement 3D finite element analysis
- Improve computational efficiency for large-scale problems

## References
See `Report.pdf` for a detailed explanation of the theory, implementation, and validation steps.

