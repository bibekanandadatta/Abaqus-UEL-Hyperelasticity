# Abaqus UEL Hyperelasticity

 
This repository contains the Fortran source code for finite strain elasticity (hyperelasticity) user element (UEL) subroutine and example input files for Abaqus/Standard. Standard displacement-based element formulation has been used. The purpose of the project is to make users familiar with developing the UEL subroutine in Abaqus/Standard using a standard displacement-based element formulation for a nonlinear (large strain) solid mechanics model. This source code contains necessary subroutines related to element operations and linear algebraic calculations.



> [!WARNING]
> This repository is not meant to be a complete guideline or tutorial of Abaqus user element (UEL) subroutines. This feature is intended for advanced users, and I strongly recommend developers and users make themselves familiar with theories related to continuum mechanics and finite element analysis, Fortran programming, and Abaqus environments before using UEL. A simpler example of Abaqus UEL subroutine implementation which can also be a good starting point for new developers or users is [Abaqus-UEL-Elasticity](https://github.com/bibekananda-datta/Abaqus-UEL-Elasticity).



## Obtaining the file

If you have `git` installed, you can clone the repository to your local machine using
```bash
git clone https://github.com/bibekananda-datta/Abaqus-UEL-Hyperelasticity.git
```

You can also `fork` the repository and sync as updates are deployed, test and develop your code by creating a separate branch.

Alternatively, you can download the files in a `zip` folder in this repository using the `code` drop-down menu on the top right corner. In this approach, you will not receive any bug fixes and updates.

> [!NOTE]
> Compiling the source code requires the LAPACK library from the Intel oneMKL package. See below for the details.


## Description of the finite element formulation

Finite element formulation for large strain problems (hyperelasticity, finite viscoelasticity, etc.) can be done in the reference configuration (total Lagrangian formulation) or the current configuration (updated Lagrangian formulation). Abaqus uses a Cauchy stress-based updated Lagrangian formulation for its built-in models. In this implementation, a second Piola-Kirchhoff (PK-II) stress-based total Lagrangian formulation has been adopted. Furthermore, two different types of constitutive models for hyperelastic materials, Neo-Hookean and Arruda-Boyce, have been implemented. For both implementations, coupled strain energy density formulations were assumed with no kinematic split. On the contrary, Abaqus uses an uncoupled split of the deformation gradient and strain energy density (split into deviatoric and volumetric parts).


> [!NOTE]
> For the built-in large displacement elements, Abaqus/Standard implements an updated Lagrangian formulation with deformation gradient and strain energy density being split into deviatoric and volumetric parts.



## Description of the repository

All the source codes are located in the `src` subdirectory and the Abaqus test cases are located in the `tests` subdirectory. The documentations are available in the `docs` subdirectory. 

|   File name   |  Description  |
| ----------    | ------------  |
| `uel_nlmech_pk2.for` | is the Fortran source code that implements PK-II stress-based Total Lagrangian user element formulation for hyperelastic materials (Neo-Hookean and Arruda-Boyce). The main `UEL` subroutine was to perform all the initial checks and the calculations are performed in a subsequent subroutine. The source code includes additional subroutines with Lagrangian interpolation functions for 4 types of 2D continuum elements (Tri3, Tri6, Quad4, and Quad8) and 4 types of 3D continuum elements (Tet4, Tet10, Hex8, Hex20) and Gaussian quadratures with reduced and full integration schemes. Body force and traction boundary conditions have not been implemented in this user subroutine, however, these can be applied by overlaying standard Abaqus elements on the user elements (to be discussed in the **Visualization** section). Since Abaqus/ Viewer does not provide native support for visualizing user elements, an additional layer of elements with the same element connectivity has been created and results at the integration points of the elements are stored using the `UVARM` subroutine. |
| `<some_module>.for` | These are the utility files with different Fortran module that are included in the main source file using `include <filename.ext>` statement at the beginning of the main source code. |
| `<...>.inp` | are the example input files prepared to be executed with the user element subroutine. Since the user-defined elements share the same topology as one of the Abaqus built-in elements, those models were built in Abaqus/CAE and then exported as input files. Later those input files were modified to include keywords and data to include user element definitions, properties, and overlaying dummy elements. |
| `addElemNLMech.py` | is a Python code that modifies a simple input file and adds the overlaying dummy elements on the user elements. For complicated input files, this will not work properly and modification of this code will be required (optional). |
| `abaqus_v6.env` | is the Abaqus environment file for Windows systems which adds the additional compiling option for the Intel oneMKL library. This needs to be in the same directory as the Abaqus jobs. |
| `runAbq.ps1` | is a PowerShell file that can invoke the Abaqus solver when executed  with the user subroutine and user-specified input file from the PowerShell terminal (optional).
| `printAbq.ps1` | is a Powershell file that can print the Abaqus `.sta` file which logs the information related to the solution process (optional). |
| Hyperelastic.pdf | is the theory and algorithm documentation for the finite element formulation and constitutive models being implemented in the provided source code. |
| Abaqus Docs.pdf | is the collection of publicly available Abaqus documentation in PDF format which is related to Abaqus UEL. The web versions of these documents are available at https://help.3ds.com. |




## Modeling in Abaqus

Since the implemented user elements have the same topology as the built-in Abaqus elements, users can build a primary model in Abaqus/CAE and then export the input (`.inp`) file. Once the input file is available, as described in the Abaqus documentation, the following information needs to be modified.


### Properties

Depending on the material model, the user needs to specify two or three properties for the material as listed below.

- Shear modulus, $\mu$
- Bulk modulus, $\kappa$
- Locking stretch, $\lambda_\mathrm{L}$
- Number of integration points, `nInt`
- Material model, `matID`
- Number of post-processed variables, `nPostVars`

> [!NOTE] 
> Use `matID = 1` for the Neo-Hookean model and `matID = 2` for the Arruda-Boyce model. For the Neo-Hookean model, the locking stretch needs to be specified as zero, and For the Arruda-Boyce model, it should be a positive real number.

> [!CAUTION] 
> The standard displacement-based element formulation implemented here is known to behave poorly for quasi-incompressible cases because of volumetric locking.



### Dummy elements for visualization

To visualize the results, an additional set of built-in Abaqus elements with the same element connectivity as the user element has been created in the input file. These additional elements (so-called dummy elements) have negligible elastic properties and thus will not affect the results. If you are using a reduced integration element from the user subroutine, then use the same type of element from Abaqus as dummy elements.




## Configuring Abaqus and executing the subroutine

To run user subroutines in Abaqus, you will need to install Intel Visual Studio and Intel oneAPI package and link them with Abaqus. Follow [this blog tutorial](https://www.bibekanandadatta.com/blog/2021/link-intel-and-vs-abaqus-2020/) if you have not done it before. Additionally, [see this blog post](https://www.bibekanandadatta.com/blog/2024/lapack-Intel-Fortran-Abaqus/) to learn how to link and use the LAPACK library from the Intel oneMKL package to Abaqus user subroutines.


Make sure that the user subroutine and input file are in the same directory. Using the `Abaqus command line terminal` or `cmd terminal` or `PowerShell terminal`, you can execute the following command from the directory to execute the subroutine.

```bash
abaqus interactive double analysis ask_delete=off job=<your_job_name> input=<input_file_name.inp> user=../src/uel_nlmech_pk2.for
```
Specify the variable names (inside < >) in the above command as needed. For additional information on executing user subroutines, check the Abaqus user manual.

If you use the PowerShell-based terminal, you can also execute the subroutine by running the `runAbq.ps1` file. Make sure to check the input file name in the file.
```bash
./runAbq
```

Users can print out the Abaqus `.sta` file by executing another PowerShell file in a different terminal opened in the same directory.
```bash
./printAbq
```


## Citation

If you use this repository (documentation or source code), please consider citing this from the following:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11078801.svg)](https://doi.org/10.5281/zenodo.11078801)

APA format:
```
Datta, B. (2024, April 28). A nonlinear user element (UEL) implementation for hyperelastic materials in Abaqus/Standard. Zenodo. https://doi.org/10.5281/zenodo.11078801.
```

BibTeX:
``` bibtex
@misc{dattaNonlinearUserElement2024,
  author       = {Datta, Bibekananda},
  title        = {{A nonlinear user element (UEL) implementation for hyperelastic materials in Abaqus/Standard}},
  month        = apr,
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.11078801},
  url          = {https://doi.org/10.5281/zenodo.11078801}
}
```