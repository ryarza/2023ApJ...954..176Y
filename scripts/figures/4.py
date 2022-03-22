import numpy as np
import matplotlib.pyplot as plt
import rytools as rt
import h5py
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Input file", required = True)
args = parser.parse_args()

rt.plot.plottex()
cgs = rt.units.get_cgs()

# Planet masses (cgs)
mps = np.linspace(10, 100, 200) * cgs['MJUP']
r_of_m = rt.planets.mass_radius_relation()
rps = r_of_m(mps / cgs['MJUP']) * cgs['RJUP']

d = rt.h5_to_dict(args.inputdata)

# Take only profiles in which stellar radius increases
mask = np.array([], dtype = int)
mask = np.append(mask, 0)
for i in range(1, len(d['max_stellar_radii'])):
    if d['max_stellar_radii'][i] > d['max_stellar_radii'][mask[-1]]:
        mask = np.append(mask, i)

d = rt.apply_mask_to_dict(d, mask)

assert ( np.diff(d['max_stellar_radii']) > 0 ).all(), np.where(np.diff(d['max_stellar_radii']) < 0)[0]

x, y, z = rt.plot.pointsToPcolormeshArgs(d['max_stellar_radii'], mps / cgs['MJUP'], d['final_rp_over_ra_sqs'])
fig, ax = plt.subplots()
im = ax.pcolormesh(x, y, np.log10(z), edgecolors = 'face', vmin = 0, vmax = 2, rasterized = True)
cbar = plt.colorbar(im)
cbar.set_label(r'$2\log_{10}\lp R_\text{SB} / R_a \rp$', rotation=90)
cbar.ax.tick_params(axis='y', direction='in')

unique, unique_idxs = np.unique(d['minimum_ejection_mass_idxs'], return_index = True)
ejec_masses = mps[d['minimum_ejection_mass_idxs']]
mass_filter = ejec_masses[unique_idxs] < 100 * cgs['MJUP']
ax.plot(d['max_stellar_radii'][unique_idxs][mass_filter], ejec_masses[unique_idxs][mass_filter] / cgs['MJUP'], color = 'black', ls = 'dashed')
ax.text(70, 85, 'Can eject envelope', ha='center', fontsize = 'x-large')
ax.text(25, 30, "Can't eject envelope", ha='center', fontsize = 'x-large')
ax.set_xscale('log')
ax.set_xlim(left = 10)
ax2 = ax.twiny()
ax2.tick_params(which = 'both', direction = 'in')
ax2.set_xlim(10 * cgs['RSUN'] / cgs['AU'], ax.get_xlim()[1] * cgs['RSUN'] / cgs['AU'])
ax2.set_xscale('log')
ax2.set_xlabel(r'$a$ $\ls\unit{\astronomicalunit}\rs$', labelpad = 10)

ax.set_ylim(10, 100)
ax.tick_params(which = 'both', direction = 'in')
ax.set_xlabel(r'$R_\star / R_\odot$')
ax.set_ylabel(r'$M_\text{SB} / M_\text{Jup}$')
plt.tight_layout()
plt.savefig('4.pdf', bbox_inches = 'tight', dpi = 300)
plt.close()
