import numpy as np
import matplotlib.pyplot as plt
import argparse
import rytools as rt
import matplotlib

rt.plot.plottex()
cgs = rt.units.get_cgs(source = 'MESA')

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Directory with data", required = True)
args = parser.parse_args()

masses = np.array([1, 10, 50, 80], dtype = int)
r_of_m = rt.planets.mass_radius_relation()
radii = r_of_m(masses) * cgs['RJUP']
models = ['total']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
lss = ['dashed', 'solid']

data = {}
data['analytical'] = {}
data['numerical'] = {}

data['analytical'] = np.empty(len(masses), dtype = dict)
data['numerical'] = np.empty(len(masses), dtype = dict)
for i in range(len(masses)):
    data['analytical'][i] = {}
    data['numerical'][i] = {}

nsamples = [-1, -1, -1, -1]
samples = [False] * 4

# Load data and compute
for midx, mass in enumerate(masses):
    for model_idx, model in enumerate(models):
        sim_path = args.inputdata + "inspirals/analytical-drag/%s/" % mass
        data['analytical'][midx][model] = rt.planets.load_orbint_output(sim_path + 'output-%s.txt' % model, sim_path + 'scalars.txt', profile_path = args.inputdata + '10rsun_profile.data', sample = samples[midx], nsamples = nsamples[midx], compute_derived = {'deltae_from_forces': False, 'L_emergent': False})
        sim_path = args.inputdata + "inspirals/numerical-drag/%s/" % mass
        data['numerical'][midx][model] = rt.planets.load_orbint_output(sim_path + 'output-%s.txt' % model, sim_path + 'scalars.txt', profile_path = args.inputdata + '10rsun_profile.data', sample = samples[midx], nsamples = nsamples[midx], compute_derived = {'deltae_from_forces': False, 'L_emergent': False})

# Plot
fig, ax = plt.subplots()
for midx, mass in enumerate(masses):
    for model_idx, model in enumerate(models):
        if model_idx == len(models) - 1: legend = str(mass) + r'$M_\text{Jup}$'
        else: legend = None
        d = data['analytical'][midx][model]
        ax.plot(d['r'][d['survived']] / cgs['RSUN'], - d['edot_from_forces']['total'][d['survived']], color = colors[midx], ls = 'dashed')
        ax.plot(d['r'][d['destroyed']['total']] / cgs['RSUN'], - d['edot_from_forces']['total'][d['destroyed']['total']], color = colors[midx], alpha = 0.2, ls = 'dashed')
        d = data['numerical'][midx][model]
        ax.plot(d['r'][d['survived']] / cgs['RSUN'], - d['edot_from_forces']['total'][d['survived']], color = colors[midx], label = legend)
        ax.plot(d['r'][d['destroyed']['total']] / cgs['RSUN'], - d['edot_from_forces']['total'][d['destroyed']['total']], color = colors[midx], alpha = 0.2)


ax.set_xscale('log')

ax.set_xlim(1e-1, 1e1)
ax.set_ylim(1e36, 1e44)

handles, labels = ax.get_legend_handles_labels()
line = matplotlib.lines.Line2D([0], [0], color='black', label = 'Analytical drag')
handles.insert(0,line) 
line = matplotlib.lines.Line2D([0], [0], color='black', ls = 'dashed', label = 'Numerical drag')
handles.insert(0,line) 

ax.legend(handles = handles, loc = 'lower left', ncol = 2)
ax.tick_params(which = 'both', direction = 'in')
ax.set_yscale('log')
ax.set_xlabel(r'$a/R_\odot$')
ax.set_ylabel(r'$\dot{E}$ $\ls\si{\erg\per\second}\rs$')
plt.tight_layout()
plt.savefig('9.pdf', bbox_inches = 'tight')
plt.close()
