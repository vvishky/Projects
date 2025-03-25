### Arguments
1. **FilePrefixStr**: Prefix of the input file names (e.g., `Gr122`).
2. **MatPtStr**: Material point identifier string (e.g., `MatPt`).
3. **FileExtension**: File extension of the input files (e.g., `dat`).
4. **NumMaterialPoints**: Number of material points (an integer, e.g., `10`).

### Example
```bash
python MohrCircle.py Gr10 MatPt dat 1
```
This command processes 10 files (e.g., `Gr122_MatPt1.dat` to `Gr122_MatPt10.dat`), generates Mohr's Circle plots for each material point, and creates additional plots for simple and volume-weighted averaged stress tensors.

## Input File Format
Each input file contains:
- First line: A comment line with the volume of the material point (e.g., `# AtomicVolume: 1.5`).
- Lines 2-4: A 3x3 stress tensor in MPa (space-separated floats), e.g.:
```
# AtomicVolume: 1.5
100.0 20.0 30.0
20.0 150.0 40.0
30.0 40.0 200.0
```

## Functionality
The script includes the following required functions:

1. **`readData(fname)`**  
   - **Input**: Filename (`fname`) of a data file.
   - **Output**: Stress tensor (`Sigma`) as a NumPy array and volume (`V`) as a float.
   - Reads the volume from the first line and the 3x3 stress tensor from subsequent lines.

2. **`AvgStress(n, SigArr, VolumeArr)`**  
   - **Inputs**: 
     - `n`: Number of material points to average (general, not hardcoded).
     - `SigArr`: Array of stress tensors.
     - `VolumeArr`: Array of corresponding volumes.
   - **Output**: Averaged stress tensor.
   - Computes a simple average (equal volumes) or volume-weighted average:
     \[
     \bar{\sigma} = \frac{\sum_{i=1}^n \sigma_i V_i}{\sum_{i=1}^n V_i}
     \]
   - Supports averaging over subsets of data (n < total points).

3. **`EigVal(Sigma)`**  
   - **Input**: Stress tensor (`Sigma`).
   - **Output**: Sorted eigenvalues of the stress tensor.
   - Uses NumPy's `linalg.eigvals` to compute eigenvalues.

4. **`plotMohrsCircle(SigmaEigval, foutNamePNG)`**  
   - **Inputs**: 
     - `SigmaEigval`: Eigenvalues of a stress tensor.
     - `foutNamePNG`: Output PNG filename.
   - Computes centers and radii of three Mohr's Circles:
     - \( C_1 = \frac{\sigma_2 + \sigma_3}{2}, R_1 = \frac{\sigma_2 - \sigma_3}{2} \)
     - \( C_2 = \frac{\sigma_1 + \sigma_3}{2}, R_2 = \frac{\sigma_1 - \sigma_3}{2} \)
     - \( C_3 = \frac{\sigma_1 + \sigma_2}{2}, R_3 = \frac{\sigma_1 - \sigma_2}{2} \)
   - Converts stresses from MPa to GPa (divides by 1000).
   - Plots circles using Matplotlib's `Circle` patch, with one shaded region.
   - Labels axes as \(\sigma\) [GPa] and \(\tau\) [GPa].
   - Saves the plot as a PNG file.

5. **`main(FilePrefixStr, MatPtStr, FileExtension, NumMaterialPoints)`**  
   - Orchestrates the workflow:
     - Reads data for all material points.
     - Computes simple average (equal volumes) and volume-weighted average stress tensors.
     - Generates Mohr's Circle plots for each material point and both averages.

## Output
The script generates PNG files:
- Individual material points: `MohrCirc_<FilePrefixStr>_<MatPtStr><n>.PNG` (e.g., `MohrCirc_Gr122_MatPt1.PNG`).
- Simple average: `MohrCirc_<FilePrefixStr>_Avg.PNG` (e.g., `MohrCirc_Gr122_Avg.PNG`).
- Volume-weighted average: `MohrCirc_<FilePrefixStr>_VolAvg.PNG` (e.g., `MohrCirc_Gr122_VolAvg.PNG`).

Each plot shows three Mohr's Circles with normal stress (\(\sigma\)) on the x-axis and shear stress (\(\tau\)) on the y-axis, both in GPa.

## Notes
- Stress values are scaled from MPa (input) to GPa for plotting.
- The script is generic and works with any number of material points, not hardcoded to 10.
- Code is well-commented to explain logic, avoiding trivial comments.
- Mohr's Circles are plotted with one filled circle (largest radius) and two outlined circles for clarity.

## Execution Environment
- Ensure a Python virtual environment with required libraries is activated.
- The script must be a standalone `.py` file, executable from the command line (not a Jupyter notebook).
