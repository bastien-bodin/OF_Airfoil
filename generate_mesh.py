import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import gmsh
import sys


def get_thickness_distribution(x: np.ndarray, max_thickness: float) -> np.ndarray:
    """
    Calculates the half-thickness distribution (y_t) for a NACA 4-digit airfoil.
    
    Uses standard NACA formula coefficients to generate a thickness distribution
    that provides the characteristic aerodynamic shape of the airfoil.
    
    Args:
        x: Normalized x-coordinates along the chord (0 to 1)
        max_thickness: Maximum thickness of the airfoil as fraction of chord (0 to 1)
    
    Returns:
        Half-thickness distribution perpendicular to the chord
    """
    coeff_glob = 5.0 * max_thickness
    term1 =  0.2969 * np.sqrt(x)
    term2 = -0.1260 * x
    term3 = -0.3516 * x**2.0
    term4 =  0.2843 * x**3.0
    term5 = -0.1015 * x**4.0
    
    y_thickness = coeff_glob * (term1 + term2 + term3 + term4 + term5)
    return y_thickness


def get_camber_line_and_gradient(x: np.ndarray, max_camber: float, camber_pos: float):
    """
    Calculates the mean camber line (y_c) and its analytical gradient (dy_c/dx).
    
    The camber line defines the curvature of the airfoil. It is calculated in two sections:
    - Forward section (0 <= x <= p): parabola passing through the origin
    - Aft section (p < x <= 1): parabola reaching the trailing edge
    
    Args:
        x: Normalized x-coordinates along the chord (0 to 1)
        max_camber: Maximum camber as fraction of chord (m)
        camber_pos: Position of maximum camber as fraction of chord (p)
    
    Returns:
        Tuple (y_camber, dyc_dx):
            - y_camber: Ordinates of the camber line
            - dyc_dx: Gradient of the camber line (analytical derivative)
    """
    y_camber = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    if max_camber == 0.0 or camber_pos == 0.0:
        return y_camber, dyc_dx

    mask_fore = (x >= 0.0) & (x <= camber_pos)
    mask_aft = (x > camber_pos) & (x <= 1.0)

    y_camber[mask_fore] = (max_camber / camber_pos**2) * (
        2.0 * camber_pos * x[mask_fore] - x[mask_fore]**2
    )
    dyc_dx[mask_fore] = (2.0 * max_camber / camber_pos**2) * (
        camber_pos - x[mask_fore]
    )

    y_camber[mask_aft] = (max_camber / (1.0 - camber_pos)**2) * (
        (1.0 - 2.0 * camber_pos) + 2.0 * camber_pos * x[mask_aft] - x[mask_aft]**2
    )
    dyc_dx[mask_aft] = (2.0 * max_camber / (1.0 - camber_pos)**2) * (
        camber_pos - x[mask_aft]
    )

    return y_camber, dyc_dx


def get_naca_coords(x: np.ndarray, digit_string: str):
    """
    Generates upper and lower surface coordinates for a NACA 4-digit airfoil.
    
    Process:
    1. Parse NACA string (e.g., "2412" -> m=0.02, p=0.4, t=0.12)
    2. Calculate camber line and its gradient
    3. Calculate thickness distribution
    4. Project thickness perpendicular to the camber line
    
    Args:
        x: Normalized x-coordinates along the chord (0 to 1)
        digit_string: 4-digit NACA code (e.g., "2412", "0012")
            - 1st digit: maximum camber in % of chord
            - 2nd digit: position of max camber in tenths of chord
            - Last 2 digits: maximum thickness in % of chord
    
    Returns:
        Tuple (x_upper, y_upper, x_lower, y_lower):
            - x_upper, y_upper: Upper surface (extrados) coordinates
            - x_lower, y_lower: Lower surface (intrados) coordinates
    """
    max_camber = float(digit_string[0]) / 100.0
    camber_pos = float(digit_string[1]) / 10.0
    max_thickness = float(digit_string[2:]) / 100.0

    y_camber, dyc_dx = get_camber_line_and_gradient(x, max_camber, camber_pos)
    y_thickness = get_thickness_distribution(x, max_thickness)
    theta = np.arctan(dyc_dx)

    x_upper = x - y_thickness * np.sin(theta)
    y_upper = y_camber + y_thickness * np.cos(theta)

    x_lower = x + y_thickness * np.sin(theta)
    y_lower = y_camber - y_thickness * np.cos(theta)
    
    return x_upper, y_upper, x_lower, y_lower


def visualize_mesh_matplotlib():
    """
    Extracts and visualizes the 2D mesh loaded in the Gmsh API with Matplotlib.
    
    This function must be called BEFORE gmsh.finalize(). It retrieves all 2D
    elements (triangles, quadrangles) from the mesh and displays them as
    polygons with their edges.
    
    The mesh is displayed with:
    - White faces and black edges
    - Preserved aspect ratio (no distortion)
    - Automatic zoom on the domain
    """
    print("Extracting mesh for visualization...")
    
    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    coords = np.array(nodeCoords).reshape(-1, 3)
    tag2idx = {tag: i for i, tag in enumerate(nodeTags)}
    
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim=2)
    
    polygons = []
    
    for i, eType in enumerate(elemTypes):
        tags = elemNodeTags[i]
        num_nodes_per_elem = gmsh.model.mesh.getElementProperties(eType)[3]
        tags_reshaped = np.array(tags).reshape(-1, num_nodes_per_elem)
        
        for row_tags in tags_reshaped:
            poly_coords = [coords[tag2idx[tag]][:2] for tag in row_tags]
            polygons.append(poly_coords)

    fig, ax = plt.subplots(figsize=(10, 8))
    
    collection = PolyCollection(polygons, edgecolors='black', facecolors='white', linewidths=0.5, alpha=0.9)
    ax.add_collection(collection)
    
    all_x = coords[:, 0]
    all_y = coords[:, 1]
    
    ax.set_xlim(all_x.min(), all_x.max())
    ax.set_ylim(all_y.min(), all_y.max())
    ax.set_aspect('equal')
    ax.set_title("Mesh Visualization (Gmsh -> Matplotlib)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    
    ax.set_xlim(-0.1, 1.1); ax.set_ylim(-0.5, 0.5)
    
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.show()


def get_cosine_spacing(n_points: int) -> np.ndarray:
    """
    Generates x-coordinates with clustering at both ends (0 and 1).
    
    Uses a cosine transformation to concentrate points at the leading edge
    and trailing edge, where gradients are highest. This improves mesh
    quality and accuracy of CFD calculations.
    
    Args:
        n_points: Number of points to generate
    
    Returns:
        Array of x-coordinates spaced according to cosine distribution (0 to 1)
    """
    beta = np.linspace(0, np.pi, n_points)
    x_distribution = 0.5 * (1.0 - np.cos(beta))
    return x_distribution


def generate_openfoam_mesh(xu, yu, xl, yl, filename="airfoil.msh"):
    """
    Generates a 3D structured C-grid mesh for OpenFOAM simulation.
    
    This function creates a mesh suitable for flow simulation around a NACA
    airfoil using Gmsh. The domain is a C-Grid type with:
    - Refined zone around the airfoil
    - Coarse zone in the far field
    - 2D -> 3D extrusion for OpenFOAM compatibility
    
    Physical groups created (for OpenFOAM boundary conditions):
    - internalField: computational volume
    - inlet: flow inlet (arcs + sides)
    - outlet: flow outlet (downstream face)
    - airfoil: airfoil wall (no-slip condition)
    - frontAndBack: front/back faces (empty type for 2D)
    
    Args:
        xu, yu: (x, y) coordinates of the upper surface (extrados)
        xl, yl: (x, y) coordinates of the lower surface (intrados)
        filename: Output mesh file name (.msh)
    
    Note:
        Mesh parameters (Lc_airfoil, R_far, L_wake) are hardcoded
        but can be modified in the function body.
    """
    gmsh.initialize()
    gmsh.model.add("naca_airfoil")

    Lc_airfoil = 0.005
    Lc_farfield = 0.5
    z_thick = 0.1
    
    R_far = 20.0
    L_wake = 30.0

    pts_upper = []
    pts_lower = []

    le_tag = gmsh.model.geo.addPoint(xu[0], yu[0], 0, Lc_airfoil)
    
    pts_upper.append(le_tag)
    pts_lower.append(le_tag)

    for i in range(1, len(xu)):
        tag = gmsh.model.geo.addPoint(xu[i], yu[i], 0, Lc_airfoil)
        pts_upper.append(tag)
    
    for i in range(1, len(xl)):
        tag = gmsh.model.geo.addPoint(xl[i], yl[i], 0, Lc_airfoil)
        pts_lower.append(tag)

    c_upper = gmsh.model.geo.addSpline(pts_upper)
    c_lower = gmsh.model.geo.addSpline(pts_lower)

    p_center = pts_upper[0]
    
    p_far_top = gmsh.model.geo.addPoint(0, R_far, 0, Lc_farfield)
    p_far_bot = gmsh.model.geo.addPoint(0, -R_far, 0, Lc_farfield)
    p_far_mid = gmsh.model.geo.addPoint(-R_far, 0, 0, Lc_farfield)
    
    p_wake_top = gmsh.model.geo.addPoint(L_wake, R_far, 0, Lc_farfield)
    p_wake_bot = gmsh.model.geo.addPoint(L_wake, -R_far, 0, Lc_farfield)

    arc_top = gmsh.model.geo.addCircleArc(p_far_mid, p_center, p_far_top)
    arc_bot = gmsh.model.geo.addCircleArc(p_far_mid, p_center, p_far_bot)
    
    l_top = gmsh.model.geo.addLine(p_far_top, p_wake_top)
    l_bot = gmsh.model.geo.addLine(p_far_bot, p_wake_bot)
    l_outlet = gmsh.model.geo.addLine(p_wake_top, p_wake_bot)

    te_upper_tag = pts_upper[-1]
    te_lower_tag = pts_lower[-1]
    l_te = gmsh.model.geo.addLine(te_upper_tag, te_lower_tag)

    loop_far = gmsh.model.geo.addCurveLoop([l_top, l_outlet, -l_bot, -arc_bot, arc_top])
    loop_airfoil = gmsh.model.geo.addCurveLoop([c_upper, l_te, -c_lower])

    surf_2d = gmsh.model.geo.addPlaneSurface([loop_far, loop_airfoil])

    new_entities = gmsh.model.geo.extrude([(2, surf_2d)], 0, 0, z_thick, [1], [1])
    
    gmsh.model.geo.synchronize()

    vol_tag = new_entities[1][1]
    surf_top_tag = new_entities[0][1]
    surf_bot_tag = surf_2d
    
    f_far_top = new_entities[2][1]
    f_outlet  = new_entities[3][1]
    f_far_bot = new_entities[4][1]
    f_arc_bot = new_entities[5][1]
    f_arc_top = new_entities[6][1]

    idx_airfoil_start = 2 + 5 
    f_airfoil_upper = new_entities[idx_airfoil_start][1]
    f_te            = new_entities[idx_airfoil_start + 1][1]
    f_airfoil_lower = new_entities[idx_airfoil_start + 2][1]

    gmsh.model.addPhysicalGroup(3, [vol_tag], name="internalField")
    
    inlet_faces = [f_far_top, f_far_bot, f_arc_bot, f_arc_top]
    gmsh.model.addPhysicalGroup(2, inlet_faces, name="inlet")
    
    gmsh.model.addPhysicalGroup(2, [f_outlet], name="outlet")
    
    wall_faces = [f_airfoil_upper, f_airfoil_lower, f_te]
    gmsh.model.addPhysicalGroup(2, wall_faces, name="airfoil")
    
    gmsh.model.addPhysicalGroup(2, [surf_top_tag, surf_bot_tag], name="frontAndBack")

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(filename)
    
    print(f"Mesh generated successfully: {filename}")

    visualize_mesh_matplotlib()
    gmsh.finalize()


if __name__ == "__main__":
    x = get_cosine_spacing(101)
    x_upper, y_upper, x_lower, y_lower = get_naca_coords(x, "2412")
    
    generate_openfoam_mesh(x_upper, y_upper, x_lower, y_lower)
    
    plt.plot(x_upper, y_upper, 'x-')
    plt.plot(x_lower, y_lower, 'x-')
    plt.axis("equal")
    plt.show()