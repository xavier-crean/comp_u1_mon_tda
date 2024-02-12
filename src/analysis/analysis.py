# # Betti Number Analysis
# ## of Monopole Current Networks in Compact $U(1)$ Lattice Gauge Theory 
# ---

# Import general libraries and set options
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import matplotlib

# Create ../../reports/
try:
    os.makedirs("../../reports/")
except FileExistsError:
    # directory already exists
    pass

# Create ../../reports/figures/
try:
    os.makedirs("../../reports/figures/")
except FileExistsError:
    # directory already exists
    pass

# Set plotting options
plt.rcParams['font.size'] = '10'
plt.rcParams["figure.figsize"] = (5,5)
cmap = matplotlib.colormaps['viridis']
cm_subsection = np.array(np.linspace(0, 255, 7), dtype=int)
colors = [ cmap.colors[x] for x in cm_subsection ]
fmt = ['_','*','p','D','^','s','o']

# Lattice Sizes
Ls = np.array([6,7,8,9,10,11,12])

# Number of Samples
N = 200

# Set number of bootstraps and random seed
N_bootstraps = 500
bootstrap_seed = 223

# Beta values for each lattice size
betas6 = np.array([0.9   , 0.95  , 0.97  , 0.975 , 0.98  , 0.99  , 0.991 , 0.992 ,
       0.993 , 0.994 , 0.995 , 0.996 , 0.997 , 0.998 , 0.999 , 1.    ,
       1.001 , 1.0016, 1.0017, 1.0018, 1.002 , 1.003 , 1.004 , 1.005 ,
       1.006 , 1.007 , 1.008 , 1.009 , 1.01  , 1.011 , 1.012 , 1.013 ,
       1.014 , 1.015 , 1.02  , 1.025 , 1.03  , 1.05  , 1.1   ])
betas7 = np.array([0.95  , 0.97  , 0.975 , 0.98  , 0.99  , 0.995 , 0.996 , 0.997 ,
       0.998 , 0.999 , 1.    , 1.001 , 1.002 , 1.0025, 1.003 , 1.0035,
       1.0037, 1.004 , 1.0045, 1.005 , 1.0057, 1.0058, 1.0059, 1.006 ,
       1.0065, 1.007 , 1.008 , 1.009 , 1.01  , 1.011 , 1.012 , 1.013 ,
       1.014 , 1.015 , 1.02  , 1.025 , 1.03  , 1.05  , 1.1   ])
betas8 = np.array([0.95  , 0.97  , 0.975 , 0.98  , 0.985 , 0.99  , 0.995 , 0.996 ,
       0.997 , 0.998 , 0.999 , 1.    , 1.001 , 1.002 , 1.003 , 1.004 ,
       1.005 , 1.006 , 1.0065, 1.007 , 1.0073, 1.0074, 1.0075, 1.008 ,
       1.0083, 1.0085, 1.009 , 1.0092, 1.0094, 1.01  , 1.011 , 1.012 ,
       1.013 , 1.014 , 1.015 , 1.02  , 1.025 , 1.03  , 1.05  ])
betas9 = np.array([0.985 , 0.99  , 0.995 , 1.    , 1.005 , 1.006 , 1.0068, 1.007 ,
       1.0072, 1.0074, 1.0076, 1.0078, 1.008 , 1.0082, 1.0084, 1.0086,
       1.0088, 1.009 , 1.0092, 1.0093, 1.0094, 1.0095, 1.0096, 1.01  ,
       1.0102, 1.0104, 1.0106, 1.0108, 1.011 , 1.0112, 1.0114, 1.0116,
       1.0118, 1.012 , 1.015 , 1.02  , 1.025 , 1.03  , 1.035 ])
betas10 = np.array([0.985 , 0.99  , 0.995 , 1.    , 1.005 , 1.006 , 1.007 , 1.0072,
       1.0074, 1.0076, 1.0078, 1.008 , 1.0082, 1.0084, 1.0086, 1.0088,
       1.009 , 1.0092, 1.0093, 1.0094, 1.0095, 1.0096, 1.0098, 1.01  ,
       1.0102, 1.0104, 1.0106, 1.0108, 1.011 , 1.0112, 1.0114, 1.0116,
       1.0118, 1.012 , 1.015 , 1.02  , 1.025 , 1.03  , 1.035 ])
betas11 = np.array([0.996 , 0.998 , 1.    , 1.002 , 1.004 , 1.006 , 1.007 , 1.0075,
        1.008 , 1.0085, 1.0088, 1.009 , 1.0091, 1.0092, 1.0093, 1.0094,
        1.0095, 1.0096, 1.0097, 1.0098, 1.0099, 1.01  , 1.0101, 1.0102,
        1.0103, 1.0104, 1.0106, 1.0108, 1.011 , 1.0115, 1.012 , 1.0125,
        1.013 , 1.014 , 1.016 , 1.018 , 1.02  , 1.022 , 1.024 ])
betas12 = np.array([0.996 , 0.998 , 1.    , 1.002 , 1.004 , 1.006 , 1.007 , 1.008 ,
       1.0082, 1.0084, 1.0086, 1.0088, 1.009 , 1.0092, 1.0094, 1.0096,
       1.0098, 1.0099, 1.01  , 1.0101, 1.0102, 1.0103, 1.0104, 1.0105,
       1.0106, 1.0108, 1.011 , 1.0112, 1.0114, 1.0116, 1.0118, 1.012 ,
       1.013 , 1.014 , 1.016 , 1.018 , 1.02  , 1.022 , 1.024])
betas = [betas6,betas7,betas8,betas9,betas10,betas11,betas12]

# ---

# ## Scatter Plot of Betti Number Observables with Error Bars
# (Error bars are computed via bootstrapping method)

# Load/Compute Betti number data

# Functions for computing Betti numbers from persistence diagram format from giotto-tda
def multiplicity(diagram,homology_dimensions=None):
    '''
    Computes multiplicity of unique points in a persistence diagram
    '''
    diagram = diagram[diagram[:, 0] != diagram[:, 1]]
    if homology_dimensions is None:
        homology_dimensions = np.unique(diagram[:, 2])
    c = []
    for dim in homology_dimensions:
        subdiagram = diagram[diagram[:, 2] == dim]
        unique, inverse, counts = np.unique(
                subdiagram, axis=0, return_inverse=True, return_counts=True
                )
        c.append(counts)
    return c

def multiplicity_total(pd):
    '''
    Computes 0-th and 1st Betti numbers
    (Note in each persistence diagram there is one finite-death time point with high multiplicity and infinite-death time points with multiplicity of the 4-torus;
      in this trivial filtration the betti numbers are simply the sum of the multiplicities)
    '''
    temp = [multiplicity(pd[s], homology_dimensions=[0,1]) for s in range(N)]
    return np.array([[np.sum(temp[s][i]) for i in range(2)] for s in range(N)])

# Load Betti number .h5 file if exists (else compute and save)
filename = f"../../data/intermediate/betti_num_{Ls[0]}-{Ls[-1]}.h5"
if os.path.exists(filename):
    with h5py.File(filename, 'r') as hf:
        h = hf['h'][()]
else:
    # Load persistence diagram .h5 files, compute Betti numbers and save as .h5 file
    pd = [[h5py.File(f'../../data/observables/pd/{Ls[l]}.{Ls[l]}.{Ls[l]}.{Ls[l]}/pers_mon_Ns={Ls[l]}{Ls[l]}{Ls[l]}{Ls[l]}_b={b:.4f}.h5', "r").get('persistence')[()] 
            for b in betas[l]] 
        for l in range(len(Ls))]
    h = np.array([[multiplicity_total(pd[l][b]) 
                   for b in range(len(betas[l]))] 
                  for l in range(len(Ls))])
    with h5py.File(filename, 'w') as hf:
        hf.create_dataset("h", data=h)


# Error bars are computed via the bootstrapping method

# Generate bootstraps
np.random.seed(bootstrap_seed)
bootstraps = np.array([[[np.sort(np.random.choice(range(N), size=N, replace=True)) for i in range(N_bootstraps)]
                    for b in range(len(betas[l]))]
                   for l in range(len(Ls))])


# Plot density of connected components $\rho_{b_{0}}$

# Generate b0 bootstrap distribution
b0_bootstrap = np.array([[np.array([h[l,:,:,0][b,bootstraps[l][b][i]] for b in range(len(betas[l]))])
                        for l in range(len(Ls))]
                       for i in range(N_bootstraps)])

plt.figure()
for l in range(len(Ls)):
    x = betas[l]
    y = np.mean(h[l,:,:,0]/(Ls[l]**4), axis=1)
    plt.scatter(x, y, s=15, marker=fmt[l], label=r'$L = '+str(int(Ls[l]))+'$', color=colors[l])
    yerr = np.std(np.mean(b0_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    plt.errorbar(x, y, ls = "None", yerr=yerr, capsize=2, color=colors[l])
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\langle \rho_{b_{0}} \rangle$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.locator_params(axis='x', nbins=10)
plt.legend(loc='upper left')
plt.savefig("../../reports/figures/p_b0_scatter.png")
plt.clf()

plt.figure()
for l in range(len(Ls)):
    x = betas[l]
    y = np.mean(h[l,:,:,0]/(Ls[l]**4), axis=1)
    plt.scatter(x, y, s=15, marker=fmt[l], label=r'$L = '+str(int(Ls[l]))+'$', color=colors[l])
    yerr = np.std(np.mean(b0_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    plt.errorbar(x, y, ls = "None", yerr=yerr, capsize=2, color=colors[l])
plt.xlim([0.99,1.02])
plt.ylim([0.007,0.019])
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\langle \rho_{b_{0}} \rangle$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.locator_params(axis='x', nbins=10)
plt.legend(loc='upper left')
plt.savefig("../../reports/figures/p_b0_scatter_ZOOM.png")
plt.clf()

# Generate b1 bootstrap distribution
b1_bootstrap = np.array([[np.array([h[l,:,:,1][b,bootstraps[l][b][i]] for b in range(len(betas[l]))])
                        for l in range(len(Ls))]
                       for i in range(N_bootstraps)])

plt.figure()
for l in range(len(Ls)):
    x = betas[l]
    y = np.mean(h[l,:,:,1]/(Ls[l]**4), axis=1)
    plt.scatter(x, y, s=15, marker=fmt[l], label=r'$L = '+str(int(Ls[l]))+'$', color=colors[l])
    yerr = np.std(np.mean(b1_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    plt.errorbar(x, y, ls = "None", yerr=yerr, capsize=4, color=colors[l])
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\langle \rho_{b_{1}} \rangle$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.locator_params(axis='x', nbins=10)
plt.legend(loc='upper right')
plt.savefig("../../reports/figures/p_b1_scatter.png")
plt.clf()

plt.figure()
for l in range(len(Ls)):
    x = betas[l]
    y = np.mean(h[l,:,:,1]/(Ls[l]**4), axis=1)
    plt.scatter(x, y, s=15, marker=fmt[l], label=r'$L = '+str(int(Ls[l]))+'$', color=colors[l])
    yerr = np.std(np.mean(b1_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    plt.errorbar(x, y, ls = "None", yerr=yerr, capsize=2, color=colors[l])
plt.xlim([0.99,1.02])
plt.ylim([0.022,0.058])
plt.xlabel(r'$\beta$')
plt.ylabel(r'$\langle \rho_{b_{1}} \rangle$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.locator_params(axis='x', nbins=10)
plt.legend(loc='upper right')
plt.savefig("../../reports/figures/p_b1_scatter_ZOOM.png")
plt.clf()

# Plot $\rho_{b_{0}}$ vs $\rho_{b_{1}}$
plt.figure()
for l in range(len(Ls)):
    x = np.mean((h[l,:,:,0]/(Ls[l]**4)) , axis=1)
    y = np.mean((h[l,:,:,1]/(Ls[l]**4)) , axis=1)
    plt.scatter(x, y, s=10, color=colors[l], marker=fmt[l], label=r'$L = '+str(int(Ls[l]))+'$')
    xerr = np.std(np.mean(b0_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    yerr = np.std(np.mean(b1_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    plt.errorbar(x, y, ls = "None", xerr=xerr, yerr=yerr, capsize=3, color=colors[l])
plt.xlabel(r'$ \langle \rho_{b_{0}} \rangle $')
plt.ylabel(r'$ \langle \rho_{b_{1}} \rangle $')
plt.tight_layout()
plt.locator_params(axis='x', nbins=10)
plt.legend(loc='upper right')
plt.savefig("../../reports/figures/p_b0_vs_p_b1_scatter.png")
plt.clf()

plt.figure()
for l in range(len(Ls)):
    x = np.mean((h[l,:,:,0]/(Ls[l]**4)) , axis=1)
    y = np.mean((h[l,:,:,1]/(Ls[l]**4)) , axis=1)
    plt.scatter(x, y, s=5, color=colors[l], marker=fmt[l], label=r'$L = '+str(int(Ls[l]))+'$')
    xerr = np.std(np.mean(b0_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    yerr = np.std(np.mean(b1_bootstrap[:,l]/(Ls[l]**4), axis=2), axis=0)
    plt.errorbar(x, y, ls = "None", xerr=xerr, yerr=yerr, capsize=2, color=colors[l])   
plt.xlabel(r'$ \langle \rho_{b_{0}} \rangle $')
plt.ylabel(r'$ \langle \rho_{b_{1}} \rangle $')
plt.xlim([0.005,0.0275])
plt.ylim([0.005,0.0275])
plt.locator_params(axis='x', nbins=5)
plt.locator_params(axis='y', nbins=5)
plt.tight_layout()
plt.legend(loc='upper left')
plt.savefig("../../reports/figures/p_b0_vs_p_b1_scatter_ZOOM.png")
plt.clf()

# ---

# ## Compute Pseudo-critical $\beta_{c}(L)$s

# For 
# 1. Average action observable $E$
# 2. Density of zeroth Betti number $\rho_{b_{0}}$
# 3. Density of first Betti number $\rho_{b_{1}}$

# Load action .h5 files and compute $\Delta S$

# Load action .h5 files
s = np.array([[h5py.File(f'../../data/observables/action/{Ls[l]}.{Ls[l]}.{Ls[l]}.{Ls[l]}/action_mon_Ns={Ls[l]}{Ls[l]}{Ls[l]}{Ls[l]}_b={betas[l][b]:.4f}.h5', "r").get('action')[()] 
      for b in range(len(betas[l]))] 
     for l in range(len(Ls))])

# Compute delta_S for use in multiple histogram reweighting
s_bar = [np.mean(s[l], axis=0) for l in range(len(Ls))]
delta_s = [s[l] - s_bar[l][np.newaxis, :] for l in range(len(Ls))]

# Import pymbar library (and turn off unnecessary warnings)
import pymbar
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias.*')


def reweight(O,E,rw_betas,l):
    '''
    Outputs reweighted mean and variance of observable O
    Uses the pymbar library which implements mulitple histogram reweighting
    '''
    rw_means = []
    rw_vars = []
    potentials = np.outer(-1*betas[l], E)
    mbar = pymbar.MBAR(potentials, np.array([N for b in betas[l]]), initialize='BAR')
    potentials = np.outer(-1*rw_betas, E)
    A = np.concatenate(O)
    mean = mbar.compute_expectations(A, u_kn=potentials, compute_uncertainty=False)['mu']
    mean2 = mbar.compute_expectations(A**2, u_kn=potentials)['mu']
    var = mean2 - mean**2        
    rw_means.append(mean)
    rw_vars.append(var)
    return rw_means, rw_vars


# Set reweighting windows and partition
n_part = 100
rw_betas_linear6 = np.linspace(0.997, 1.006, n_part)
rw_betas_linear7 = np.linspace(1.003, 1.008, n_part)
rw_betas_linear8 = np.linspace(1.005 , 1.009, n_part)
rw_betas_linear9 = np.linspace(1.007, 1.0104, n_part)
rw_betas_linear10 = np.linspace(1.008, 1.0104, n_part)
rw_betas_linear11 = np.linspace(1.009, 1.0104, n_part)
rw_betas_linear12 = np.linspace(1.0096, 1.0104, n_part)
rw_betas_linear = [rw_betas_linear6,rw_betas_linear7,rw_betas_linear8,rw_betas_linear9,rw_betas_linear10,rw_betas_linear11,rw_betas_linear12]


# Reweight $E$

# For each lattice size, load/compute multiple histogram reweighting of b0 observable
rw_means_linear = []
rw_vars_linear = []
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_linear_rw{np.min(rw_betas_linear[l])}-{np.max(rw_betas_linear[l])}_Ls={Ls[l]}_E.h5"
    if os.path.exists(filename):
        with h5py.File(filename, 'r') as hf:
            rw_means_linear.append(hf['rw_means_linear'][()])
            rw_vars_linear.append(hf['rw_vars_linear'][()])
    else:
        rw_means_linear_temp, rw_vars_linear_temp = reweight(s[l]/(6*(Ls[l]**4)), delta_s[l], rw_betas_linear[l], l)
        with h5py.File(filename, 'w') as hf:
            hf.create_dataset("rw_means_linear", data=rw_means_linear_temp)
            hf.create_dataset("rw_vars_linear", data=rw_vars_linear_temp)
        rw_means_linear.append(rw_means_linear_temp)
        rw_vars_linear.append(rw_vars_linear_temp)
rw_means_linear_E = np.concatenate(rw_means_linear)
rw_vars_linear_E = np.concatenate(rw_vars_linear)


# Reweight $\rho_{b_{0}}$

# For each lattice size, load/compute multiple histogram reweighting of b0 observable
rw_means_linear = []
rw_vars_linear = []
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_linear_rw{np.min(rw_betas_linear[l])}-{np.max(rw_betas_linear[l])}_Ls={Ls[l]}_H0.h5"
    if os.path.exists(filename):
        with h5py.File(filename, 'r') as hf:
            rw_means_linear.append(hf['rw_means_linear'][()])
            rw_vars_linear.append(hf['rw_vars_linear'][()])
    else:
        rw_means_linear_temp, rw_vars_linear_temp = reweight(h[l,:,:,0]/(Ls[l]**4), delta_s[l], rw_betas_linear[l], l)
        with h5py.File(filename, 'w') as hf:
            hf.create_dataset("rw_means_linear", data=rw_means_linear_temp)
            hf.create_dataset("rw_vars_linear", data=rw_vars_linear_temp)
        rw_means_linear.append(rw_means_linear_temp)
        rw_vars_linear.append(rw_vars_linear_temp)
rw_means_linear_b0 = np.concatenate(rw_means_linear)
rw_vars_linear_b0 = np.concatenate(rw_vars_linear)


# Reweight $\rho_{b_{1}}$

# For each lattice size, load/compute multiple histogram reweighting of b1 observable
rw_means_linear = []
rw_vars_linear = []
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_linear_rw{np.min(rw_betas_linear[l])}-{np.max(rw_betas_linear[l])}_Ls={Ls[l]}_H1.h5"
    if os.path.exists(filename):
        with h5py.File(filename, 'r') as hf:
            rw_means_linear.append(hf['rw_means_linear'][()])
            rw_vars_linear.append(hf['rw_vars_linear'][()])
    else:
        rw_means_linear_temp, rw_vars_linear_temp = reweight(h[l,:,:,1]/(Ls[l]**4), delta_s[l], rw_betas_linear[l], l)
        with h5py.File(filename, 'w') as hf:
            hf.create_dataset("rw_means_linear", data=rw_means_linear_temp)
            hf.create_dataset("rw_vars_linear", data=rw_vars_linear_temp)
        rw_means_linear.append(rw_means_linear_temp)
        rw_vars_linear.append(rw_vars_linear_temp)
rw_means_linear_b1 = np.concatenate(rw_means_linear)
rw_vars_linear_b1 = np.concatenate(rw_vars_linear)

# Locate psuedo-critical $\beta_{c}(L)$ values as argmax of the reweighted variance curves

pseudo_bcs_E = np.array([rw_betas_linear[l][np.argmax(rw_vars_linear_E[l])]
                       for l in range(len(Ls))])
pseudo_bcs_b0 = np.array([rw_betas_linear[l][np.argmax(rw_vars_linear_b0[l])]
                       for l in range(len(Ls))])
pseudo_bcs_b1 = np.array([rw_betas_linear[l][np.argmax(rw_vars_linear_b1[l])]
                       for l in range(len(Ls))])


# #### Compute error bars

# Generate Delta S bootstrap distribution
delta_s_bootstrap = np.array([[np.array([delta_s[l][b][bootstraps[l][b][i]]
                           for b in range(len(betas[l]))])
                          for l in range(len(Ls))]
                         for i in range(N_bootstraps)])

# Generate E bootstrap distribution
E_bootstrap = np.array([[np.array([(s[l]/(6*(Ls[l]**4)))[b,bootstraps[l][b][i]] for b in range(len(betas[l]))])
                        for l in range(len(Ls))]
                       for i in range(N_bootstraps)])


# For $E$, compute multiple histogram reweighting of bootstraps

# For each lattice size, compute and/or multiple histogram reweighting of E bootstraps
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_data_bootstraps_seed={bootstrap_seed}_rw_Ls={Ls[l]}_E.h5"
    if not os.path.exists(filename):
        rw_means_linear_E_bs_l = []
        rw_vars_linear_E_bs_l = []
        for i in range(N_bootstraps):
            rwm, rwv = reweight(E_bootstrap[i][l], delta_s_bootstrap[i][l], rw_betas_linear[l],l)
            rw_means_linear_E_bs_l.append(rwm)
            rw_vars_linear_E_bs_l.append(rwv)
        with h5py.File(filename, 'w') as hf:
            hf.create_dataset("rw_means_linear_E_bs", data=rw_means_linear_E_bs_l)
            hf.create_dataset("rw_vars_linear_E_bs", data=rw_vars_linear_E_bs_l)
rw_means_linear_E_bs = []
rw_vars_linear_E_bs = []
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_data_bootstraps_seed={bootstrap_seed}_rw_Ls={Ls[l]}_E.h5"
    with h5py.File(filename, 'r') as hf:
        rw_means_linear_E_bs.append(hf['rw_means_linear_E_bs'][()])
        rw_vars_linear_E_bs.append(hf['rw_vars_linear_E_bs'][()])
rw_means_linear_E_bs = np.concatenate(rw_means_linear_E_bs,axis=1)
rw_vars_linear_E_bs = np.concatenate(rw_vars_linear_E_bs,axis=1)


# For $\rho_{b_{0}}$, compute multiple histogram reweighting of bootstraps

# For each lattice size, compute and/or multiple histogram reweighting of b0 bootstraps
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_data_bootstraps_seed={bootstrap_seed}_rw_Ls={Ls[l]}_b0.h5"
    if not os.path.exists(filename):
        rw_means_linear_b0_bs_l = []
        rw_vars_linear_b0_bs_l = []
        for i in range(N_bootstraps):
            rwm, rwv = reweight(b0_bootstrap[i][l], delta_s_bootstrap[i][l], rw_betas_linear[l],l)
            rw_means_linear_b0_bs_l.append(rwm)
            rw_vars_linear_b0_bs_l.append(rwv)
        with h5py.File(filename, 'w') as hf:
            hf.create_dataset("rw_means_linear_b0_bs", data=rw_means_linear_b0_bs_l)
            hf.create_dataset("rw_vars_linear_b0_bs", data=rw_vars_linear_b0_bs_l)
rw_means_linear_b0_bs = []
rw_vars_linear_b0_bs = []
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_data_bootstraps_seed={bootstrap_seed}_rw_Ls={Ls[l]}_b0.h5"
    with h5py.File(filename, 'r') as hf:
        rw_means_linear_b0_bs.append(hf['rw_means_linear_b0_bs'][()])
        rw_vars_linear_b0_bs.append(hf['rw_vars_linear_b0_bs'][()])
rw_means_linear_b0_bs = np.concatenate(rw_means_linear_b0_bs,axis=1)
rw_vars_linear_b0_bs = np.concatenate(rw_vars_linear_b0_bs,axis=1)


# For $\rho_{b_{1}}$, compute multiple histogram reweighting of bootstraps


# For each lattice size, compute and/or multiple histogram reweighting of b1 bootstraps
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_data_bootstraps_seed={bootstrap_seed}_rw_Ls={Ls[l]}_b1.h5"
    if not os.path.exists(filename):
        rw_means_linear_b1_bs_l = []
        rw_vars_linear_b1_bs_l = []
        for i in range(N_bootstraps):
            rwm, rwv = reweight(b1_bootstrap[i][l], delta_s_bootstrap[i][l], rw_betas_linear[l],l)
            rw_means_linear_b1_bs_l.append(rwm)
            rw_vars_linear_b1_bs_l.append(rwv)
        with h5py.File(filename, 'w') as hf:
            hf.create_dataset("rw_means_linear_b1_bs", data=rw_means_linear_b1_bs_l)
            hf.create_dataset("rw_vars_linear_b1_bs", data=rw_vars_linear_b1_bs_l)
rw_means_linear_b1_bs = []
rw_vars_linear_b1_bs = []
for l in range(len(Ls)):
    filename = f"../../data/intermediate/intermediate_data_bootstraps_seed={bootstrap_seed}_rw_Ls={Ls[l]}_b1.h5"
    with h5py.File(filename, 'r') as hf:
        rw_means_linear_b1_bs.append(hf['rw_means_linear_b1_bs'][()])
        rw_vars_linear_b1_bs.append(hf['rw_vars_linear_b1_bs'][()])
rw_means_linear_b1_bs = np.concatenate(rw_means_linear_b1_bs,axis=1)
rw_vars_linear_b1_bs = np.concatenate(rw_vars_linear_b1_bs,axis=1)


pseudo_bcs_E_bs = np.array([[rw_betas_linear[l][np.argmax(rw_vars_linear_E_bs[i][l])]
                             for l in range(len(Ls))]
                            for i in range(N_bootstraps)])
yerr_E = np.sqrt(np.var(pseudo_bcs_E_bs, axis=0))

pseudo_bcs_b0_bs = np.array([[rw_betas_linear[l][np.argmax(rw_vars_linear_b0_bs[i][l])]
                             for l in range(len(Ls))]
                            for i in range(N_bootstraps)])
yerr_b0 = np.sqrt(np.var(pseudo_bcs_b0_bs, axis=0))

pseudo_bcs_b1_bs = np.array([[rw_betas_linear[l][np.argmax(rw_vars_linear_b1_bs[i][l])]
                             for l in range(len(Ls))]
                            for i in range(N_bootstraps)])
yerr_b1 = np.sqrt(np.var(pseudo_bcs_b1_bs, axis=0))


# ### Print pseudo-critical betas

# Save as a .csv file


# Save .csv withOUT header or first column
csv = np.array(
    [[pseudo_bcs_E[l],yerr_E[l],pseudo_bcs_b0[l],yerr_b0[l],pseudo_bcs_b1[l],yerr_b1[l]] for l in range(len(Ls))]
)
np.savetxt("../../reports/pseudo_bcs.csv", csv, delimiter=",",fmt='%s')


# Save .csv with header and first column
head = np.array(['L','E','E_err','rho_b0','rho_b0_err','rho_b1','rho_b1_err']).reshape(-1,7)
csv = np.array(
    [[f'{Ls[l]}',f'{pseudo_bcs_E[l]}',f'{yerr_E[l]}',f'{pseudo_bcs_b0[l]}',f'{yerr_b0[l]}',f'{pseudo_bcs_b1[l]}',f'{yerr_b1[l]}'] for l in range(len(Ls))]
)
csv_head = np.concatenate([head,csv],axis=0)
np.savetxt("../../reports/pseudo_bcs_HEADER.csv", csv_head, delimiter=",", fmt='%s')


# Save as a .tex file

def format_NIST(m,u):
    '''Formats measurement and standard error in NIST:
    https://physics.nist.gov/cgi-bin/cuu/Info/Constants/definitions.html
    '''
    m = float(m)
    u = float(u)
    u_order = int(np.floor(np.log10(u)))
    round_m = '{:.{}f}'.format(m,-1*u_order)
    round_u = str(int(np.round(u * 10**(-u_order))))
    return round_m + '(' + round_u + ')'


# Write table with pseudo_bcs to .tex file 
filename = '../../reports/pseudo_beta_c.tex'
with open(filename,'w') as f:
    data = [(format_NIST(pseudo_bcs_E[l],yerr_E[l]),format_NIST(pseudo_bcs_b0[l],yerr_b0[l]), format_NIST(pseudo_bcs_b1[l],yerr_b1[l])) for l in range(len(Ls))]
    tab_content = r'''
        $6$ & ${}$ & ${}$ & ${}$ \\
        $7$ & ${}$ & ${}$ & ${}$ \\
        $8$ & ${}$ & ${}$ & ${}$ \\
        $9$ & ${}$ & ${}$ & ${}$ \\
        $10$ & ${}$ & ${}$ & ${}$ \\
        $11$ & ${}$ & ${}$ & ${}$ \\
        $12$ & ${}$ & ${}$ & ${}$
    '''.format(*[item for sublist in data for item in sublist])
    tab_content.replace('\n','')
    f.write(tab_content)


# # Finite-size scaling analysis

# Perform a finite-size scaling analysis according to the asymptotic relation $\beta_{c}(L) = \beta_{c} + \sum_{k=1}^{\infty} B_{k} L^{-4k}$ where sum truncation is set by $k_{max} < \infty$


# Import libraries from scikit-learn for polynomial regression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression


# Set spacetime dimension
d = 4
inverted_d = 1/d

# Set x as inverse power of the volume V = L^d
x = Ls**-(1/inverted_d)

def ffs_rchi2(Ls_for_reg,deg,y,y_bs):
    '''Function to output range,k_max,reduced_chi2,beta_c,beta_c_err'''
    # Perform polynomial regression
    poly_features = PolynomialFeatures(degree=deg, include_bias=False)
    X_poly = poly_features.fit_transform(x.reshape(-1,1))
    lin_reg = LinearRegression().fit(X_poly[Ls_for_reg], y[Ls_for_reg])
    
    # Compute intercept error and coefficient error
    intercepts = [LinearRegression().fit(X_poly[Ls_for_reg], y_bs[i][Ls_for_reg]).intercept_ 
                   for i in range(N_bootstraps)]
    intercept_err = np.sqrt(np.var(intercepts))
    coefs = [LinearRegression().fit(X_poly[Ls_for_reg], y_bs[i][Ls_for_reg]).coef_ 
                       for i in range(N_bootstraps)]
    coef_err = np.sqrt(np.var(coefs))
    
    # Compute reduced chi-squared statistic
    f_obs = y[Ls_for_reg]
    f_exp = lin_reg.predict(X_poly[Ls_for_reg])
    var = np.var(y_bs[:,Ls_for_reg], axis=0)
    chi2 = np.sum((f_obs-f_exp)**2 / var )
    nu = f_obs.size - (deg+1)
    reduced_chi2 = chi2 / nu
    
    return f'{Ls[Ls_for_reg]}',deg,f'{reduced_chi2:.3f}',lin_reg.intercept_,intercept_err,format_NIST(lin_reg.intercept_,intercept_err)


# ### Print FSS results for $E$

# Save as a .csv file


# Save .csv with header or first column
range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
ls = [ffs_rchi2(range_[i],d,pseudo_bcs_E,pseudo_bcs_E_bs)[:-1] for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
csv = np.array([element for tup in ls for element in tup]).reshape(-1,5)
np.savetxt("../../reports/ffs_E.csv", csv, delimiter=",",fmt='%s')


# Save .csv with header
head = np.array(['Range','k_max','chi_per_dof','beta_c','beta_c_err'])
range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
ls = [ffs_rchi2(range_[i],d,pseudo_bcs_E,pseudo_bcs_E_bs)[:-1] for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
csv = np.concatenate([np.array(['Range','k_max','chi_per_dof','beta_c','beta_c_err']).reshape(-1,5),np.array([element for tup in ls for element in tup]).reshape(-1,5)],axis=0)
np.savetxt("../../reports/ffs_E_HEADER.csv", csv, delimiter=",",fmt='%s')


# Save as .tex file

# Write table with pseudo_bcs to .tex file 
filename = '../../reports/ffs_E.tex'
with open(filename,'w') as f:
    range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
    ls = [ffs_rchi2(range_[i],d,pseudo_bcs_E,pseudo_bcs_E_bs) for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
    data = np.array([element for tup in ls for element in tup]).reshape(-1,6)[:,[1,2,5]]
    tab_content = r'''
        $10,11,12$ & ${}$ & ${}$ & ${}$ \\
        \hline
        $8,9,10,11,12$ & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
        \hline
        $6,7,8,9,10,11,12$ & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$
    '''.format(*[item for sublist in data for item in sublist])
    tab_content.replace('\n','')
    f.write(tab_content)


# ### Plot Polynomial Regression for $\rho_{b_{0}}$

# ### Print FSS results for $\rho_{b_{0}}$

# Save as a .csv file


# Save .csv withOUT header
range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
ls = [ffs_rchi2(range_[i],d,pseudo_bcs_b0,pseudo_bcs_b0_bs)[:-1] for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
csv = np.array([element for tup in ls for element in tup]).reshape(-1,5)
np.savetxt("../../reports/ffs_rho_b0.csv", csv, delimiter=",",fmt='%s')


# Save .csv with header
head = np.array(['Range','k_max','chi_per_dof','beta_c','beta_c_err'])
range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
ls = [ffs_rchi2(range_[i],d,pseudo_bcs_b0,pseudo_bcs_b0_bs)[:-1] for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
csv = np.concatenate([np.array(['Range','k_max','chi_per_dof','beta_c','beta_c_err']).reshape(-1,5),np.array([element for tup in ls for element in tup]).reshape(-1,5)],axis=0)
np.savetxt("../../reports/ffs_rho_b0_HEADER.csv", csv, delimiter=",",fmt='%s')


# Save as .tex file


# Write table with pseudo_bcs to .tex file 
filename = '../../reports/ffs_rho_b0.tex'
with open(filename,'w') as f:
    range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
    ls = [ffs_rchi2(range_[i],d,pseudo_bcs_b0,pseudo_bcs_b0_bs) for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
    data = np.array([element for tup in ls for element in tup]).reshape(-1,6)[:,[1,2,5]]
    tab_content = r'''
        $10,11,12$ & ${}$ & ${}$ & ${}$ \\
        \hline
        $8,9,10,11,12$ & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
        \hline
        $6,7,8,9,10,11,12$ & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$
    '''.format(*[item for sublist in data for item in sublist])
    tab_content.replace('\n','')
    f.write(tab_content)


# ### Plot Polynomial Regression for $\rho_{b_{1}}$

# ### Print FSS results for $\rho_{b_{1}}$

# Save as a .csv file


# Save .csv withOUT header
range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
ls = [ffs_rchi2(range_[i],d,pseudo_bcs_b1,pseudo_bcs_b1_bs)[:-1] for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
csv = np.array([element for tup in ls for element in tup]).reshape(-1,5)
np.savetxt("../../reports/ffs_rho_b1.csv", csv, delimiter=",",fmt='%s')


# Save .csv with header
head = np.array(['Range','k_max','chi_per_dof','beta_c','beta_c_err'])
range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
ls = [ffs_rchi2(range_[i],d,pseudo_bcs_b1,pseudo_bcs_b1_bs)[:-1] for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
csv = np.concatenate([np.array(['Range','k_max','chi_per_dof','beta_c','beta_c_err']).reshape(-1,5),np.array([element for tup in ls for element in tup]).reshape(-1,5)],axis=0)
np.savetxt("../../reports/ffs_rho_b1_HEADER.csv", csv, delimiter=",",fmt='%s')


# Save as .tex file

# Write table with pseudo_bcs to .tex file 
filename = '../../reports/ffs_rho_b1.tex'
with open(filename,'w') as f:
    range_ = [[4,5,6],[2,3,4,5,6],[0,1,2,3,4,5,6]]
    ls = [ffs_rchi2(range_[i],d,pseudo_bcs_b1,pseudo_bcs_b1_bs) for i in range(len(range_)) for d in range(1,len(range_[i])-1)]
    data = np.array([element for tup in ls for element in tup]).reshape(-1,6)[:,[1,2,5]]
    tab_content = r'''
        $10,11,12$ & ${}$ & ${}$ & ${}$ \\
        \hline
        $8,9,10,11,12$ & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
        \hline
        $6,7,8,9,10,11,12$ & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$ \\
               & ${}$ & ${}$ & ${}$
    '''.format(*[item for sublist in data for item in sublist])
    tab_content.replace('\n','')
    f.write(tab_content)


