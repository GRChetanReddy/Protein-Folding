import matplotlib.pyplot as plt
import numpy as np 
from zimm_bragg import ZimmBragg, ZimmBragg_2s
from zimm_bragg import plot_Gprofile



#########################################################
################## ZIMM-BRAGG ###########################
#########################################################

NRES = 15
SIGMA = 1e-3

############### SCENARIO A S >> 1 ########################

s_high = np.linspace(1, 2, 11)
t_values = np.linspace(278, 368, 6)

# compare s values for T = 278 K
G_highs_T278 = np.zeros((len(s_high), NRES+1))
for i, s in enumerate(s_high):
    zb_a1 = ZimmBragg(s, SIGMA, NRES)
    G_highs_T278[i, :] = zb_a1.free_energy_profile(temps=[278])
plot_Gprofile(G_highs_T278, output="ZB_A_varyS", title="Scenario A: vary s", 
              labels=[f"{s:.1f}" for s in s_high], nres=NRES)

# compare temperatures at s = 1.5
zb_a2 = ZimmBragg(1.5, SIGMA, NRES)
G_T_s10 = zb_a2.free_energy_profile(temps=t_values)
plot_Gprofile(G_T_s10, title="Scenario A: vary T", output="ZB_A_varyT", 
              labels=[f"{t:.0f} K" for t in t_values], nres=NRES)

##########################################################

############### SCENARIO B S << 1 ########################

s_low = np.linspace(0.1, 1.0, 10)
t_values = np.linspace(278, 368, 6)

# compare s values for T = 278 K
G_lows_T278 = np.zeros((len(s_low), NRES+1))
for i, s in enumerate(s_low):
    zb_b1 = ZimmBragg(s, SIGMA, NRES)
    G_lows_T278[i, :] = zb_b1.free_energy_profile(temps=[278])
plot_Gprofile(G_lows_T278, output="ZB_B_varyS", title="Scenario B: vary s", 
              labels=[f"{s:.1f}" for s in s_low])

# compare temperatures at s = 0.1
zb_b2 = ZimmBragg(0.1, SIGMA, NRES)
G_T_s01 = zb_b2.free_energy_profile(temps=t_values)
plot_Gprofile(G_T_s01, title="Scenario B: vary T", output="ZB_B_varyT", 
              labels=[f"{t:.0f} K" for t in t_values])

##########################################################

############## SCENARIO C S_h > 1, S_c < 1 ###############
s_h = 1.5
s_c = 0.1
t_values = np.linspace(278, 368, 6)
# 1 for helix, 0 for helix
residues = np.array([0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0])

zb_c1 = ZimmBragg_2s(s_h, s_c, residues, SIGMA, NRES)
G_2s_T = zb_c1.free_energy_profile(temps=t_values)  
plot_Gprofile(G_2s_T, title=f"Scenario C: s={s_h}, {s_c}", output="ZB_C",
              labels=[f"{t:.0f} K" for t in t_values])

##########################################################

################# Fraction Helicity ######################

t_values = np.linspace(100, 2000, 50)
s_values = np.logspace(-5, 5, 11)
s_lab = np.linspace(-5, 5, 11)
print("s     Melting Temperatures")
for s, sl in zip(s_values, s_lab):
    zb = ZimmBragg(s, SIGMA)
    fracs, melt_temp = zb.thermal_unfolding_curves(t_values)
    plt.plot(t_values, fracs, label = f"10$^{{{sl:.0f}}}$")
    print(s, melt_temp)
fracs, melt_temp = zb_c1.thermal_unfolding_curves(t_values)
print("2S", melt_temp)
plt.plot(t_values, fracs, label="2S model")
plt.legend()
plt.title("Thermal Unfolding Curves")
plt.xlabel("T [K]")
plt.ylabel(r"$\theta$")
plt.savefig("./plots/ThermalUnfolding.png")
plt.show()
