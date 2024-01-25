import numpy as np

"""
parameter part
"""

PHI_0 = 2.067833848 * 1e-15
k_B = 1.38e-23
T = 0.5
k_BT = k_B * T



I_p_1, I_p_2 = 2e-6 , 2e-6  # Amp
I_m_1, I_m_2 = 7e-9, 7e-9                                # Amp
R_1, R_2 = 371, 371                                # ohm
C_1, C_2 = 4e-9, 4e-9                              # F
L_1, L_2 = 1e-9, 1e-9                              # H

quick_doubler = lambda x1, x2: np.hstack([np.array([x1] * 2), np.array([x2]*2)])
I_p, I_m = quick_doubler(I_p_1, I_p_2), quick_doubler(I_m_1, I_m_2)

m_c = C_1
m_1 = C_1
m_2 = C_2
x_c = PHI_0 / (2 * np.pi)
time_scale_factor = 1
t_c = np.sqrt(L_1 * C_1)


U0_1 = m_c * x_c**2 / t_c
U0_2 = m_2 * x_c**2 / t_c
kappa_1, kappa_2, kappa_3, kappa_4 = k_BT/U0_1, k_BT/U0_1, k_BT/U0_1, k_BT/U0_1




lambda_1 = 2 * np.sqrt(L_1 * C_1) / (C_1 * R_1)
theta_1  = 1
eta_1    = (L_1 * C_1 / R_1**2)**0.25 * np.sqrt(2 * kappa_1)

lambda_2 = 2 * np.sqrt(L_1 * C_1) / (C_2 * R_2)
theta_2  = 1 / (C_2/C_1)
eta_2    = (L_1 * C_1 / R_1**2)**0.25 * np.sqrt(2 * kappa_2 * (R_1 * C_1) / (R_2 * C_2))

lambda_3 = 2 * np.sqrt(L_1 * C_1) / (C_1 * R_1)
theta_3  = 4
eta_3    = (L_1 * C_1 / R_1**2)**0.25 * np.sqrt(2 * kappa_3)

lambda_4 = 2 * np.sqrt(L_1 * C_1) / (C_2 * R_2)
theta_4  = 4 / (C_2/C_1)
eta_4    = (L_1 * C_1 / R_1**2)**0.25 * np.sqrt(2 * kappa_4 * (R_1 * C_1) / (R_2 * C_2))

gamma = 20


beta_1 = 2 * np.pi * L_1 * I_p_1 / PHI_0; 
beta_2 = 2 * np.pi * L_2 * I_p_2 / PHI_0;

d_beta_1 = 2 * np.pi * L_1 * I_m_1 / PHI_0; 
d_beta_2 = 2 * np.pi * L_2 * I_m_2 / PHI_0;


_lambda = np.array([lambda_1, lambda_2, lambda_3, lambda_4])
_theta  = np.array([theta_1, theta_2, theta_3, theta_4])
_eta  =   np.array([eta_1, eta_2, eta_3, eta_4])