# 2D Finite Element Stress Analysis

This software is designed for conducting stress analysis in a 2D plane using the plane strain approximation method. It offers a comprehensive set of tools to investigate stress distribution in materials, which is invaluable in engineering and materials science.

## Features

- **Isoparametric Shape Functions**: The software utilizes shape and gradient defining functions to establish the shape and gradient information for stress analysis.
- **Mesh Formation**: Enables the creation of a mesh with defined parameters (number of elements, nodes, and dimensions) to represent the structure.
- **Material Model - Plane Strain**: Utilizes material properties such as Young's Modulus and Poisson's ratio to model materials under plane strain conditions.
- **Global Stiffness Matrix Calculation**: Computes the stiffness matrix using numerical integration for each element based on its connectivity and node coordinates.
- **Nodal Forces and Boundary Conditions**: Assigns nodal forces and boundary conditions, crucial for simulating real-world scenarios.
- **Solving Linear System**: Solves the linear system to obtain the displacement values.
- **Displacement Visualization**: Offers visualization of the displacement using contour plots, aiding in understanding stress distribution.

## Usage

To utilize the software:
1. Ensure necessary libraries (NumPy, Matplotlib) are installed.
2. Define mesh and material properties based on your analysis requirements.
3. Run the script to generate the stress analysis results.

## Output
![Strain visualisation](https://raw.githubusercontent.com/YoussefNassar-1959/2d_FEM_Solver/main/Stress%20visualization.png)

## Installation

To get started with the software:

1. Clone or download this repository.
2. Ensure you have the required libraries installed (`numpy`, `matplotlib`).
3. Execute the script following the usage instructions mentioned in the repository.

## Contributing

Contributions are encouraged! Whether you want to report issues, suggest improvements, or submit pull requests, your input is valued. Please follow these guidelines:
- Open a clear and descriptive issue for bug reports or enhancement suggestions.
- Fork the repository, make changes, and submit a pull request with concise details of the modifications.

## License

This software is under [XYZ License](license-file), permitting users to [describe the permissions and restrictions].

## Acknowledgments

Acknowledgments to contributors, libraries, or any resources used in developing this software.
