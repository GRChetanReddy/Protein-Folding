import numpy as np
from scipy.special import binom
import matplotlib.pyplot as plt
import itertools
import os

'''
Assumptions: 
    free energies for nucleaation and coil-helix transition are not temperature dependent
    T dependency through weights, s and sigma
    partition sum is T dependent
'''

class ZimmBragg:
    
    def __init__(self, s, sigma, nres=15, ref_temp = 298):
        self.R = 8.314463    # J/(K mol)
        # G of h->c and nucleation (unit: J/mol)
        self.G_hc = -self.R*ref_temp*np.log(s)
        self.G_nuc = -self.R*ref_temp*np.log(sigma)
        self.nres = nres
                
    def partition_sum(self, temp):
        s = np.exp(-self.G_hc/(self.R * temp))
        sigma = np.exp(-self.G_nuc/(self.R * temp))
        l = np.array([0, 1])
        r = np.array([[1, 1]]).T
        W = np.array([[s, 1], [sigma * s, 1]])
        W_prod = np.linalg.matrix_power(W, self.nres)
        self.Z = l@W_prod@r

    def free_energy_profile(self, temps, out=True):
        # store probabilities of macrostate nh as function of T
        self.probs = np.zeros((len(temps), self.nres+1))        
        self.G_profile = np.zeros((len(temps), self.nres+1))
        for i, temp in enumerate(temps):
            self.partition_sum(temp)
            for n in range(self.nres+1):
                G_n, prob_n = self.free_energy_nhelix(temp, n)
                self.G_profile[i, n] = G_n
                self.probs[i, n] = prob_n
        if out==True:
            return self.G_profile

    def free_energy_nhelix(self, temp, nhelix):
        if nhelix == 0:
            return 1/self.Z, 1/self.Z
        else:
            # maximum v (number of helix stretches separated by at least 1 coil)
            if self.nres % 2 == 0:
                max_v = int(self.nres/2)
            else:
                max_v = int((self.nres+1)/2)
            G_tot = 0
            prob_tot = 0
            all_gnv = binom(self.nres, nhelix)
            for v in range(1, np.min([nhelix, max_v])+1):                
                # g_nv is number of possible combinations to have nhelix in helix
                # from https://pubs.aip.org/aip/jcp/article/38/4/934/207167/On-the-Helix-Coil-Equilibrium-in-Polypeptides eq. 8, µ=1
                g_nv = binom(nhelix-1, v-1) * binom(self.nres-nhelix+1, v)
                all_gnv -= g_nv
                # G for every microstate described by (n,v)
                G_nv = nhelix * self.G_hc + v * self.G_nuc
                # statistical weight of (n,v) microstate
                prob_nv = g_nv * np.exp(-G_nv/(temp*self.R)) / self.Z
                G_tot += G_nv * prob_nv
                prob_tot += prob_nv
            if all_gnv != 0:
                print(f"Warning: Not all states considered: {all_gnv} states missing!!")
            return G_tot, prob_tot
    
    def thermal_unfolding_curves(self, temps):
        # calculate probabilities for given temperatures
        self.free_energy_profile(temps, out=False)
        # take updated probs and calculate fractional helicities
        fract_hel, melting_temp = fractional_helicities(self.probs, temps)
        return fract_hel, melting_temp
        

class ZimmBragg_2s:
    
    def __init__(self, s_h, s_c, res_class, sigma, nres=15, ref_temp=298):
        self.res_class = res_class
        self.nres = nres
        self.R = 8.314463    # J/(K mol)
        # G of h->c and nucleation (unit: J/mol), no T-dependency
        self.G_hc_h = -self.R*ref_temp*np.log(s_h)
        self.G_hc_c = -self.R*ref_temp*np.log(s_c)
        self.G_nuc = -self.R*ref_temp*np.log(sigma)
        
    def partition_sum(self, temp):
        sh = np.exp(-self.G_hc_h/(self.R * temp))
        sc = np.exp(-self.G_hc_c/(self.R * temp))
        sigma = np.exp(-self.G_nuc/(self.R * temp))
        l = np.array([0, 1])
        r = np.array([[1, 1]]).T
        W_h = np.array([[sh, 1], [sigma * sh, 1]])
        W_c = np.array([[sc, 1], [sigma * sc, 1]])
        W_prod = np.identity(2)
        for i in range(self.nres):
            if self.res_class[i] == 1:
                W_prod = W_prod@W_h
            elif self.res_class[i] == 0:
                W_prod = W_prod@W_c
            else:
                print("Warning: Only 0 and 1 allowed to classify state!")
        self.Z = l@W_prod@r
    
    def free_energy_profile(self, temps, out=True):
        
        # store probabilities of macrostate nh as function of T
        self.probs = np.zeros((len(temps), self.nres+1))
        microstates = self.microstates_data()
        
        G_profile = np.zeros((len(temps), self.nres+1))
        for i, temp in enumerate(temps):
            self.partition_sum(temp)
            for n in range(self.nres+1):
                # select the microstates which have n helices
                nh_arrays = [[microstates[:, i]] for i in range(np.shape(microstates)[1]) if microstates[0, i] == n]
                microstates_nh = np.concatenate(nh_arrays)
                G_profile[i, n], self.probs[i, n] = self.free_energy_nhelix(temp, microstates_nh)
        if out == True:
            return G_profile

    def microstates_data(self):
        # get 4 parameters of all microstates: number of helices in s_h, s_c region & 
        # number of stretches in sequence, number of helices
        all_combs = list(itertools.product([0, 1], repeat=self.nres))
        all_combs = np.array(all_combs)
        n_h = np.sum(all_combs, axis=-1)
        n_h_sh = all_combs @ self.res_class.T
        n_h_sc = n_h - n_h_sh
        all_v = np.zeros_like(all_combs[:, 0])
        for i, row in enumerate(all_combs):
            tmp = np.concatenate(([0], row, [0]))
            tmp = np.diff(tmp)
            v = np.sum(tmp == 1)
            all_v[i] = v
        params = np.array([n_h, n_h_sh, n_h_sc, all_v])
        return params

    def free_energy_nhelix(self, temp, ms_data):
        if ms_data[0, 0] == 0:
            return 1/self.Z, 1/self.Z
        else:            
            # G of every microstate
            G_params = np.array([0, self.G_hc_h, self.G_hc_c, self.G_nuc])
            G_ms = ms_data @ G_params
            #weights = self.s_h**ms_data[:, 1] * self.s_c**ms_data[:, 2] * self.sigma**ms_data[:, 3]
            prob = np.exp(-G_ms/(self.R * temp)) / self.Z
            G_tot = G_ms @ prob
            prob_tot = np.sum(prob)
            return G_tot, prob_tot    
    
    def thermal_unfolding_curves(self, temps):
        # calculate probabilities for given temperatures
        self.free_energy_profile(temps, out=False)
        # take updated probs and calculate fractional helicities
        fract_hel, melting_temp = fractional_helicities(self.probs, temps)
        return fract_hel, melting_temp

def fractional_helicities(probs, temps):
    nres_p1 = np.shape(probs)[1]
    nh = np.linspace(0, nres_p1-1, nres_p1)   
    fract_hel = probs @ nh / (nres_p1-1)
    melting_temp = "Not reached"
    for a in range(len(fract_hel)):
        if np.linalg.norm(fract_hel[a] - 0.5) < 0.1:
            melting_temp = temps[a]
            break
        
    return fract_hel, melting_temp
                
def plot_Gprofile(G_values, nres=15, title=None, labels=None, output=None):
    
    x = range(nres+1)
    if labels is None:
        plt.plot(x, G_values.T)
    else:
        for g, lab in zip(G_values, labels):
            plt.plot(x, g, label=lab)
        plt.legend()
    plt.xlabel("number of helix residues")
    plt.ylabel("G [J/mol]")
    if title:
        plt.title(title)  
    if output:
        os.makedirs('./plots', exist_ok=True)
        plt.savefig(f"./plots/{output}")
    plt.show()   
    
    
    
