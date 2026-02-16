import os


def create_numerical_settings(case_dir: str = "."):
    """
    Generates numerical schemes and solver settings for steady-state incompressible RANS.
    
    Creates two essential configuration files for simpleFoam simulations:
    
    1. fvSchemes: Discretization schemes
       - Time: steadyState (no time derivatives)
       - Convection (U): linearUpwind (2nd order, accurate for lift/drag)
       - Convection (turbulence): upwind (1st order, stable)
       - Laplacian/gradient: corrected (handles non-orthogonal Gmsh meshes)
    
    2. fvSolution: Linear solvers and SIMPLE algorithm settings
       - Pressure: GAMG (efficient multigrid solver)
       - Velocity/turbulence: smoothSolver with symGaussSeidel
       - SIMPLEC algorithm (consistent=yes) for faster convergence
       - 2 non-orthogonal correctors (important for Gmsh meshes)
    
    Args:
        case_dir: OpenFOAM case root directory (default: current directory)
    
    Note:
        These settings are optimized for external aerodynamics (airfoils).
        Relaxation factors are set to 0.9 (aggressive) due to SIMPLEC.
    """
    system_dir = os.path.join(case_dir, "system")
    
    fv_schemes_content = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{
    default         steadyState;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""

    fv_solution_content = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-6;
        relTol          0.1;
        smoother        GaussSeidel;
    }
    
    "(U|k|omega|epsilon)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
    }
}

SIMPLE
{
    nNonOrthogonalCorrectors 2;
    consistent      yes;

    residualControl
    {
        p               1e-4;
        U               1e-4;
        "(k|omega)"     1e-4;
    }
}

relaxationFactors
{
    equations
    {
        U               0.9;
        ".*"            0.9;
    }
}

// ************************************************************************* //
"""

    with open(os.path.join(system_dir, "fvSchemes"), "w") as f:
        f.write(fv_schemes_content)
        
    with open(os.path.join(system_dir, "fvSolution"), "w") as f:
        f.write(fv_solution_content)

    print(f"--> Numerical settings (fvSchemes/fvSolution) created in {system_dir}")


if __name__ == "__main__":
    create_numerical_settings()