import numpy as np
import mesa_reader as mr
import matplotlib.pyplot as plt
import rytools as rt
import argparse

rt.plot.plottex()
cgs = rt.units.get_cgs(source = 'MESA')

r_of_m = rt.planets.mass_radius_relation()

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Input file", required = True)
args = parser.parse_args()

prof = rt.stars.load_mesa_profile(args.inputdata)

# Planet masses (cgs)
mps = np.array([1, 10, 80]) * cgs['MJUP']
# Planet radii (cgs)
rps = np.array([r_of_m(mp / cgs['MJUP']) for mp in mps]) * cgs['RJUP']

core_idx = rt.stars.post_main_sequence_core_idx(rho=prof['rho'], x_h1=prof['x_h1'], criterion='hydrogen_mass_fraction')
# Binding energy from r to infinity
ebind_above = np.empty(prof['n'])
for k in range(prof['n']):
  ebind_above[k] = rt.stars.binding_energy(r = prof['r'], rho = prof['rho'], menc = prof['menc'], idx0 = k, idx1 = -1)

# Secondary axis with orbital separation
def m_of_a(a):
  try:
    if len(a) > 1:
      out = np.empty_like(a)
      for out_idx in range(len(out)):
       closest_idx, _ = rt.nearest(prof['r'], a[out_idx] * cgs['RSUN'])
       out[out_idx] = prof['menc'][closest_idx]
  except:
    out = prof['menc'][rt.nearest(prof['r'], a * cgs['RSUN'])[0]]
  return out / cgs['MSUN']

def a_of_m(mval):
  if type(mval) == np.ndarray:
    out = np.empty_like(mval)
    for out_idx in range(len(out)):
     closest_idx, _ = rt.nearest(prof['menc'], mval[out_idx] * cgs['MSUN'])
     out[out_idx] = prof['r'][closest_idx]
  else: out = prof['r'][rt.nearest(prof['menc'], mval * cgs['MSUN'])[0]]
  return out / cgs['RSUN']

fig, ax = plt.subplots()
ax.plot(prof['r'] / cgs['RSUN'], - ebind_above, label = r'$E_\text{bind}\lp>a\rp$', color = 'black')
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

for mpidx, mp in enumerate(mps):

  # Compute change in orbital energy
  eorb0 = - cgs['G'] * prof['menc'][-1] * mp / 2 / prof['r'][-1]
  eorbf = - cgs['G'] * prof['menc']     * mp / 2 / prof['r']
  delta_eorb = eorbf - eorb0
  ebind_ra = np.empty(prof['n'])

  destruction_idx_ram_pressure = rt.planets.ram_pressure_disruption_idx(prof['r'], prof['rho'], mp, rps[mpidx])
  destruction_idx_tides = rt.planets.tidal_disruption_idx(prof['r'], prof['rho'], mp, rps[mpidx])
  destruction_idx = np.amax(np.array([destruction_idx_ram_pressure, destruction_idx_tides]))

  ax.plot(prof['r'][destruction_idx+1:] / cgs['RSUN'], - delta_eorb[destruction_idx+1:], label = r'$\Delta E_\text{orb}\lp' + str(int(mp/cgs['MJUP'])) + r'M_\text{J}\rp$', color = cycle[mpidx])
  ax.plot(prof['r'][:destruction_idx+1] / cgs['RSUN'], - delta_eorb[:destruction_idx+1], ls = 'dotted', color = cycle[mpidx])

ax.set_yscale('log')
ax.set_xlim(1e-1, None)
ax.set_xscale('log')
ax.set_ylim(1e44, 1e49)
ax.set_xlabel(r'$a / R_\odot$')
ax.set_ylabel(r'$-E$ $\ls \unit{\erg} \rs$')
ax.tick_params(which = 'both', direction = 'in')

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xscale('log')
avals = np.array([1e0, 1e1, 1e2])
mvals = m_of_a(avals)
ax2.set_xticks(avals)
ax2.set_xticklabels(["%.2f" % i for i in mvals ])
ax2.tick_params(which = 'both', direction = 'in')
ax2.minorticks_off()
ax2.set_xlabel(r'$M_\text{enc} / M_\odot$', labelpad = 10)

ax.legend(loc = 0, ncol = 2)

plt.tight_layout()
plt.savefig('3.pdf', bbox_inches = 'tight')
plt.close()
