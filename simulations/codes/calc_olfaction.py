import numpy as np
import math

pN = 0.64 * 2
Z = 10
L = 2000
M = 4
tpop = 500000
delta = 1.5

nu_fmg = pN * Z / L
p_fmg_z = (1 - np.exp(-nu_fmg)) ** M
p_fmg = 1 - (1 - p_fmg_z) ** (L / Z)
print("FMG", p_fmg, p_fmg * tpop)

nu_sdg = pN * M * Z / L
p_sdg_z = 1 - np.sum(
    [((nu_sdg**k) * np.exp(-nu_sdg)) / math.factorial(k) for k in range(M)]
)
p_sdg = 1 - (1 - p_sdg_z) ** (L / Z)
print("SDG", p_sdg, p_sdg * tpop)

nu_poss = pN * delta / L
E_poss = pN * ((nu_poss) ** (M - 1))
p_poss = 1 - np.exp(-E_poss)
print("poss", p_poss, p_poss * tpop)
