import os
import pathlib
import matplotlib
import matplotlib.figure
import matplotlib.style
import numpy as np
import pandas
import unyt

os.chdir(pathlib.Path(__file__).parent.parent)
matplotlib.style.use('~/.config/matplotlib/style.mplstyle')

data = pandas.read_csv('data/planets.csv')

data = data.rename(columns={
    "pl_radj": "sb_radius",
    "pl_bmassj": "sb_mass",
    "pl_orbsmax": "orbital_separation",
    "st_mass": "star_mass"
})

data = data[data['sb_mass'] > 1]
data = data.to_dict(orient='list')

data['sb_mass'] = data['sb_mass'] * unyt.m_jup
data['sb_radius'] = data['sb_radius'] * unyt.r_jup
data['orbital_separation'] = data['orbital_separation'] * unyt.au
data['star_mass'] = data['star_mass'] * unyt.m_sun
data['keplerian_speed'] =\
    np.sqrt(unyt.G * data['star_mass'] / data['orbital_separation'])
data['sb_accretion_radius'] =\
    2 * unyt.G * data['sb_mass'] / pow(data['keplerian_speed'], 2)
data['rp_over_ra'] = data['sb_radius'] / data['sb_accretion_radius']

print(len(data['rp_over_ra'][data['rp_over_ra'] > 1]))
print(len(data['rp_over_ra'][data['rp_over_ra'] < 1]))

fig = matplotlib.figure.Figure()
axis = fig.add_subplot()
bins = np.logspace(-1, 2, 30)
axis.hist(data['rp_over_ra'], bins=bins, color='seagreen', alpha=0.8)
# ax.set_xlabel(r'$R_\text{SB}/R_a$')
axis.set_xlabel(r'Ratio between geometrical and gravitational radii')
axis.set_xscale('log')
axis.set_ylabel(r'Number of exoplanets')
axis.set_xlim(1e-1, 1e2)
axis.set_ylim(0, 100)
fig.savefig('plots/histogram.pdf')
