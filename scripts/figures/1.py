import numpy as np
import csv
import rytools
import matplotlib.pyplot as plt
import matplotlib
from collections import Counter

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Input file", required = True)
args = parser.parse_args()

rytools.plot.plottex()

cgs = rytools.units.get_cgs()

#Full file data
with open(args.inputdata) as csvfile:
    reader = csv.DictReader(csvfile, skipinitialspace=True)
    d = {name: [] for name in reader.fieldnames}
    for row in reader:
        for name in reader.fieldnames:
            d[name].append(row[name])

# New array with only masses, radius of the planet, and semimajor axis
data = {}

data['Rp']    = np.array(d['pl_radj'   ], dtype = np.float64) * cgs['RJUP']
data['Mp']    = np.array(d['pl_bmassj' ], dtype = np.float64) * cgs['MJUP']
data['a' ]    = np.array(d['pl_orbsmax'], dtype = np.float64) * cgs['AU'  ]
data['Mstar'] = np.array(d['st_mass'   ], dtype = np.float64) * cgs['MSUN']

#Derived quantities
data['vkep'] = np.sqrt( cgs['GNEWT'] * data['Mstar'] / data['a'] )
data['Ra'] = cgs['GNEWT'] * data['Mp'] / data['vkep'] / data['vkep']

#Planet mass filter
for key in ['Rp', 'a', 'Mstar', 'vkep', 'Ra']:
	data[key] = data[key][data['Mp'] >= ( 5 * cgs['MJUP'] )]
data['Mp'] = data['Mp'][data['Mp'] >= ( 5 * cgs['MJUP'] ) ]

#Histogram both a < 0.1 au and a > 0.1 au
data2 = data.copy()
for key in ['Rp', 'Mp', 'Mstar', 'vkep', 'Ra']:
	data2[key] = data[key][data['a'] <= ( 0.1 * cgs['AU'] )]
data2['a'] = data['a'][data['a'] <= ( 0.1 * cgs['AU'] ) ]

data3 = data.copy()
for key in ['Rp', 'Mp', 'Mstar', 'vkep', 'Ra']:
	data3[key] = data[key][data['a'] >= ( 0.1 * cgs['AU'] )]
data3['a'] = data['a'][data['a'] >= ( 0.1 * cgs['AU'] ) ]

fig, ax = plt.subplots()
ax.hist(2 * np.log10(data2['Rp'] / data2['Ra']), bins = np.linspace(-10, 4, 29), color = 'seagreen', alpha = 0.8, label = r'$a<\qty{e-1}{\astronomicalunit}$')
ax.hist(2 * np.log10(data3['Rp'] / data3['Ra']), bins = np.linspace(-10, 4, 29), color = 'royalblue', alpha = 0.8, label = r'$a>\qty{e-1}{\astronomicalunit}$')
#General plot settings
ax.set_xlabel(r'$2\log_{10}\left(R_\text{SB}/R_a\right)$')
ax.set_ylabel(r'$N_\text{SB}$')
ax.set_xlim(-10, 2)
ax.set_ylim(0, 22)
ax.legend(loc = 0)
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
ax.tick_params(which = 'both', direction="in")
plt.tight_layout()
plt.savefig('1.pdf', bbox_inches = 'tight')
plt.close()
