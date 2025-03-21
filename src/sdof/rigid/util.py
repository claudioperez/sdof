import numpy as np

EPS = np.finfo(float).eps

def CaySO3(w):
    theta = np.linalg.norm(w)
    if theta < 1e-10:
        return np.eye(3)
    else:
        W = HatSO3(w)
        return np.eye(3) + (np.sin(theta) / theta) * W + ((1 - np.cos(theta)) / (theta ** 2)) * (W @ W)

def ExpSO3(omega):
    """
    Exponential map from so(3) to SO(3).

    Parameters:
    omega (np.array): 3x1 vector in the Lie algebra so(3)

    Returns:
    np.array: 3x3 rotation matrix in SO(3)
    """
    theta = np.linalg.norm(omega)
    if theta < 1e-10:
        return np.eye(3)
    omega_hat = HatSO3(omega)
    return np.eye(3) + (np.sin(theta) / theta) * omega_hat + ((1 - np.cos(theta)) / (theta ** 2)) * np.dot(omega_hat, omega_hat)

def dExpSO3(phi):
    """
    Derivative of the exponential map for SO(3)
    """
    phi_norm = np.linalg.norm(phi)
    if phi_norm < 1e-10:
        return 0.5 * Spin(phi)
    else:
        axis = phi / phi_norm
        skew_phi = Spin(phi)
        a = (1 - np.cos(phi_norm)) / (phi_norm ** 2)
        b = (phi_norm - np.sin(phi_norm)) / (phi_norm ** 3)
        return (
            0.5 * skew_phi +
            a * (skew_phi @ skew_phi) +
            b * (skew_phi @ skew_phi @ skew_phi)
        )

def dExpSO3(omega):
    """
    Compute the differential of the exponential map on SO(3).

    Parameters:
    omega (np.array): 3x1 vector in the Lie algebra so(3)

    Returns:
    np.array: 3x3 matrix representing the differential of the exponential map
    """
    theta = np.linalg.norm(omega)
    if theta < 1e-10:
        return np.eye(3)

    omega_hat = HatSO3(omega)
    omega_hat_sq = omega_hat @ omega_hat

    a1 = (1 - np.cos(theta)) / (theta ** 2)
    a2 = (theta - np.sin(theta)) / (theta ** 3)

    return np.eye(3) + a1 * omega_hat + a2 * omega_hat_sq

def LogSO3(R):
    """
    Logarithm map from SO(3) to so(3) using Spurrier's algorithm.

    Parameters:
    R (np.array): 3x3 rotation matrix in SO(3)

    Returns:
    np.array: 3x1 vector in the Lie algebra so(3)
    """
    theta = np.arccos((np.trace(R) - 1) / 2)
    if theta < 1e-10:
        return np.zeros(3)
    omega_hat = (theta / (2 * np.sin(theta))) * (R - R.T)
    return VeeSO3(omega_hat)

def dLogSO3(omega):
    """
    Differential of the logarithm map (inverse of the exponential tangent) on SO(3).

    Parameters:
    R (np.array): 3x3 rotation matrix in SO(3)

    Returns:
    np.array: 3x3 matrix representing the differential of the logarithm map
    """
    # omega = LogSO3(R)
    theta = np.linalg.norm(omega)
    if theta < 1e-10:
        return np.eye(3)
    omega_hat = HatSO3(omega)
    A = (1 - np.cos(theta)) / (theta ** 2)
    B = (theta - np.sin(theta)) / (theta ** 3)
    return np.eye(3) - 0.5 * omega_hat + (1 / theta ** 2) * (1 - A / (2 * B)) * np.dot(omega_hat, omega_hat)

def HatSO3(vec):
    """
    Hat operator for so(3).

    Parameters:
    vec (np.array): 3x1 vector

    Returns:
    np.array: 3x3 skew-symmetric matrix
    """
    return np.array([
        [        0, -vec[2],  vec[1]],
        [ vec[2],         0, -vec[0]],
        [-vec[1],  vec[0],        0]
    ])

def VeeSO3(mat):
    """
    Vee operator for so(3).

    Parameters:
    mat (np.array): 3x3 skew-symmetric matrix

    Returns:
    np.array: 3x1 vector
    """
    return np.array([mat[2, 1], mat[0, 2], mat[1, 0]])
