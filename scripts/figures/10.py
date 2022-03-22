import numpy as np
import matplotlib.pyplot as plt
import rytools as rt
import argparse
import h5py

rt.plot.plottex()
cgs = rt.units.get_cgs()

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Path to input data", required = True)
args = parser.parse_args()

masses = np.array([1, 10, 50, 80], dtype = int)
r_of_m = rt.planets.mass_radius_relation()
radii = r_of_m(masses) * cgs['RJUP']
models = ['total']
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
lss = ['dashed', 'solid']


data = np.empty(len(masses), dtype = dict)
for i in range(len(data)):
    data[i] = {}

samples = [False] * 4
nsamples = [0] * 4

#compute_derived = {'deltae_from_forces': True, 'L_emergent': True}
compute_derived = False

# Load data and compute
for midx, mass in enumerate(masses):
    for model_idx, model in enumerate(models):
        sim_path = args.inputdata + "inspirals/numerical-drag/%s/" % mass
        data[midx][model] = rt.planets.load_orbint_output(sim_path + 'output-%s.txt' % model, sim_path + 'scalars.txt', profile_path = args.inputdata + '10rsun_profile.data', sample = samples[midx], nsamples = nsamples[midx], compute_derived = compute_derived)
        f = rt.h5_to_dict(sim_path + 'postprocessed.h5')
        data[midx][model]['L_emergent'] = f['L_emergent']
        data[midx][model]['L_combined'] = {}
        data[midx][model]['L_combined']['L'] = f['L_combined_L']
        data[midx][model]['L_combined']['t'] = f['L_combined_t']
#        f = h5py.File(sim_path + 'postprocessed.h5', 'w')
#        f.create_dataset('L_emergent', data = data[midx][model]['L_emergent'])
#        f.create_dataset('L_combined_L', data = data[midx][model]['L_combined']['L'])
#        f.create_dataset('L_combined_t', data = data[midx][model]['L_combined']['t'])
#        f.close()

mosaic = [['A', 'C'],
          ['B', 'C']]
fig, ax = plt.subplot_mosaic(mosaic = mosaic, figsize = (10, 4.8))

for key in ax.keys():
    ax[key].tick_params(which = 'both', direction = 'in')
    ax[key].loglog()

ax['A'].set_xlim(1e-1, 1e1)
ax['A'].set_ylim(1e-6, 1e1)
ax['B'].sharex(ax['A'])
ax['A'].tick_params('x', labelbottom=False)

for midx, mass in enumerate(masses):
    for model_idx, model in enumerate(models):
        d = data[midx][model]
        mask = d['r'] > 0
        t_insp = np.abs(d['r'] / d['rdot'] )
        ax['A'].plot(d['r'][mask] / cgs['RSUN'], t_insp[mask] / cgs['YEAR'])

        ax['B'].plot(d['r'] / cgs['RSUN'], d['L_emergent'], color = colors[midx])
        ax['B'].plot(d['r'] / cgs['RSUN'], - d['edot_from_forces']['total'], color = colors[midx], ls = 'dashed')

        ax['C'].plot(d['L_combined']['t'] / cgs['YEAR'], ( d['L_combined']['L'] + np.amax(d['prof']['L']) ) / np.amax(d['prof']['L']), label = '$' + str(mass) + r'M_\text{Jup}$')


ax['A'].plot([0], [0], label = r'$\tau_\text{insp}$', ls = 'solid', color = 'black')

ax['B'].plot([0], [0], label = r'$\Delta L$', ls = 'solid', color = 'black')
ax['B'].plot([0], [0], label = r'$\dot{E}$', ls = 'dashed', color = 'black')
 
ax['C'].set_xlabel(r'$t$ $\ls\unit{\year}\rs$')
ax['C'].set_xlim(1e-3, 1e4)
ax['C'].set_ylim(1e0, 1e5)
ax['B'].set_xlabel(r'$a/R_\odot$')
ax['A'].set_ylabel(r'$t$ $\ls\unit{\year}\rs$')
ax['B'].set_ylabel(r'$\dot{E}$ $\ls\unit{\erg\per\second}\rs$')
ax['C'].set_ylabel(r'$L/L_\star$')
ax['A'].plot(d['prof']['r'] / cgs['RSUN'], d['prof']['t_kh'] / cgs['YEAR'], ls = 'dashed', color = 'black', label = r'$\tau_\text{KH}$')

ax['A'].legend(loc = 'lower right', ncol = 2, handlelength = 1)
ax['B'].legend(loc = 0, ncol = 2, handlelength = 1)
ax['C'].legend(loc = 0, ncol = 1, handlelength = 1)

plt.tight_layout()
plt.savefig('10.pdf')
plt.close()
