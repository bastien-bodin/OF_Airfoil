import numpy as np

def get_thickness_distribution(x: np.ndarray, max_thickness: float) -> np.ndarray:
    """
    Calculates the half-thickness distribution (y_t) for a NACA 4-digit airfoil.
    """
    # Coefficients for the NACA 4-digit thickness distribution
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
    """
    y_camber = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    # Avoid division by zero for symmetric airfoils (m=0, p=0)
    if max_camber == 0.0 or camber_pos == 0.0:
        return y_camber, dyc_dx

    # Mask for the forward section (0 <= x <= p)
    mask_fore = (x >= 0.0) & (x <= camber_pos)
    # Mask for the aft section (p < x <= 1)
    mask_aft = (x > camber_pos) & (x <= 1.0)

    # Forward section calculations
    # y_c = (m / p^2) * (2px - x^2)
    # dy_c/dx = (2m / p^2) * (p - x)
    y_camber[mask_fore] = (max_camber / camber_pos**2) * (
        2.0 * camber_pos * x[mask_fore] - x[mask_fore]**2
    )
    dyc_dx[mask_fore] = (2.0 * max_camber / camber_pos**2) * (
        camber_pos - x[mask_fore]
    )

    # Aft section calculations
    # y_c = (m / (1-p)^2) * ((1-2p) + 2px - x^2)
    # dy_c/dx = (2m / (1-p)^2) * (p - x)
    y_camber[mask_aft] = (max_camber / (1.0 - camber_pos)**2) * (
        (1.0 - 2.0 * camber_pos) + 2.0 * camber_pos * x[mask_aft] - x[mask_aft]**2
    )
    dyc_dx[mask_aft] = (2.0 * max_camber / (1.0 - camber_pos)**2) * (
        camber_pos - x[mask_aft]
    )

    return y_camber, dyc_dx

def get_naca_coords(x: np.ndarray, digit_string: str):
    """
    Generates the upper and lower coordinates for a NACA 4-digit airfoil.
    """
    # Parse NACA string (e.g., "2412")
    max_camber = float(digit_string[0]) / 100.0      # m
    camber_pos = float(digit_string[1]) / 10.0       # p
    max_thickness = float(digit_string[2:]) / 100.0  # t

    # 1. Get mean camber line and its gradient
    y_camber, dyc_dx = get_camber_line_and_gradient(x, max_camber, camber_pos)

    # 2. Get thickness distribution
    y_thickness = get_thickness_distribution(x, max_thickness)

    # 3. Calculate local angle theta using the analytical gradient
    # We use arctan because we have the slope dy/dx
    theta = np.arctan(dyc_dx)

    # 4. Project thickness perpendicular to the camber line
    x_upper = x - y_thickness * np.sin(theta)
    y_upper = y_camber + y_thickness * np.cos(theta)

    x_lower = x + y_thickness * np.sin(theta)
    y_lower = y_camber - y_thickness * np.cos(theta)
    
    return x_upper, y_upper, x_lower, y_lower

import matplotlib.pyplot as plt

x = np.linspace(0,1,101)

x_upper, y_upper, x_lower, y_lower = get_naca_coords(x, "2412")

plt.plot(x_upper, y_upper)
plt.plot(x_lower, y_lower)

plt.axis("equal")

plt.show()