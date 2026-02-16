import os
import subprocess
import re


def setup_openfoam_structure(case_dir: str = "."):
    """
    Creates the necessary OpenFOAM directory structure and minimal controlDict.
    
    Sets up the standard OpenFOAM case structure with three main directories:
    - 0: Initial and boundary conditions
    - constant: Mesh and physical properties
    - system: Numerical schemes and solver settings
    
    A minimal controlDict is written to system/ with default settings for
    simpleFoam (steady-state incompressible solver).
    
    Args:
        case_dir: Root directory for the OpenFOAM case (default: current directory)
    """
    system_dir = os.path.join(case_dir, "system")
    constant_dir = os.path.join(case_dir, "constant")
    zero_dir = os.path.join(case_dir, "0")

    os.makedirs(system_dir, exist_ok=True)
    os.makedirs(constant_dir, exist_ok=True)
    os.makedirs(zero_dir, exist_ok=True)

    control_dict_content = """
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
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     simpleFoam;

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          1;

writeControl    timeStep;

writeInterval   100;

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""
    
    with open(os.path.join(system_dir, "controlDict"), "w") as f:
        f.write(control_dict_content)
    
    print(f"--> OpenFOAM structure created in {case_dir}")


def convert_msh_to_openfoam(msh_filename: str, case_dir: str = "."):
    """
    Converts a Gmsh .msh file to OpenFOAM format and patches boundary conditions for 2D.
    
    This function performs two critical steps:
    1. Runs gmshToFoam utility to convert the mesh
    2. Modifies constant/polyMesh/boundary to set correct boundary types:
       - 'frontAndBack' -> 'empty' (required for 2D simulations)
       - 'airfoil' -> 'wall' (no-slip boundary condition)
    
    Args:
        msh_filename: Path to the Gmsh mesh file (.msh)
        case_dir: OpenFOAM case directory (default: current directory)
    
    Note:
        Requires OpenFOAM environment to be sourced (gmshToFoam must be in PATH).
        The function will suppress standard output from gmshToFoam to reduce clutter.
    """
    
    print(f"--> Converting {msh_filename} to OpenFOAM format...")
    
    try:
        subprocess.run(
            ["gmshToFoam", msh_filename, "-case", case_dir], 
            check=True,
            stdout=subprocess.DEVNULL 
        )
        print("    Conversion successful.")
    except subprocess.CalledProcessError:
        print("Error: gmshToFoam failed. Make sure OpenFOAM is sourced.")
        return
    except FileNotFoundError:
        print("Error: Command 'gmshToFoam' not found.")
        return

    boundary_file = os.path.join(case_dir, "constant", "polyMesh", "boundary")
    
    if not os.path.exists(boundary_file):
        print(f"Error: {boundary_file} not found.")
        return

    print("--> Patching boundary file for 2D simulation...")
    
    with open(boundary_file, 'r') as f:
        content = f.read()

    pattern_empty = r"(frontAndBack\s*\{[^\}]*type\s+)(patch)(;)"
    
    if re.search(pattern_empty, content):
        new_content = re.sub(
            pattern_empty, 
            r"\1empty\3", 
            content
        )
        print("    'frontAndBack' set to type 'empty'.")
    else:
        print("    Warning: 'frontAndBack' patch not found or already correct.")
        new_content = content

    pattern_wall = r"(airfoil\s*\{[^\}]*type\s+)(patch)(;)"
    if re.search(pattern_wall, new_content):
         new_content = re.sub(
            pattern_wall, 
            r"\1wall\3", 
            new_content
        )
         print("    'airfoil' set to type 'wall'.")

    with open(boundary_file, 'w') as f:
        f.write(new_content)
        
    print("--> Boundary file patched successfully.")


if __name__ == "__main__":
    setup_openfoam_structure()
    
    if os.path.exists("airfoil.msh"):
        convert_msh_to_openfoam("airfoil.msh")
    else:
        print("Please generate 'airfoil.msh' first.")