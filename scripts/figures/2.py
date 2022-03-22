import numpy as np
import matplotlib.pyplot as plt
import rytools as rt
from matplotlib.ticker import FormatStrFormatter
import matplotlib

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Input file", required = True)
args = parser.parse_args()

cgs = rt.units.get_cgs(source = 'MESA')

rt.plot.plottex()

plotDir = './'
prof = rt.stars.load_mesa_profile(args.inputdata)

print("Radius [RSun]: ", prof['r'][-1] / cgs['RSUN'])
print("Luminosity [LSun]: ", prof['L'][-1] / cgs['LSUN'])

prof['vkep'] = np.sqrt( cgs['GNEWT'] * prof['menc'] / prof['r'] )
prof['mach'] = prof['vkep'] / prof['cs']

n_points = 200
min_idx, _ = rt.nearest(prof['r'], 0.1 * cgs['RSUN'])
min_idx = 0
# Range of planet masses (in Jupiter masses)
mp = np.logspace(-2, 2, n_points) * cgs['MJUP']
# Corresponding range of radii (in Jupiter radii)
r_of_m = rt.planets.mass_radius_relation()
rp = r_of_m(mp / cgs['MJUP']) * cgs['RJUP']
# vesc from the planet
vesc = np.sqrt( cgs['GNEWT'] * mp / rp )
rho_p = mp / ( 4 * np.pi * pow(rp, 3) / 3 )

jiaspruit_f = np.outer(prof['rho'] * pow(prof['vkep'], 2), pow(rho_p * pow(vesc, 2), -1 ) )
ra = 2 * cgs['GNEWT'] * np.outer(pow(prof['vkep'], -2), mp )
rp_over_ra = np.outer( prof['menc'] / prof['r'], rp / mp ) / 2
a_over_rt = np.outer(prof['r'] / pow(prof['menc'], 1/3), pow(mp, 1/3) / rp)
epsilon_rho = np.abs(np.maximum( np.outer( 1 / prof['h_rho'], rp ), ra * np.outer( 1 / prof['h_rho'], np.ones(n_points) )))

plot_x = 'a'
if plot_x == 'mass':
  plot_x_coord = prof['menc'] / cgs['MSUN']
else:
  plot_x_coord = prof['r'] / cgs['RSUN']

fields = {}
fields['eps_rho'] = {}
fields['eps_rho']['array'] = np.log10(epsilon_rho)
fields['eps_rho']['tex_name'] = r'$\log_{10}\lp\varepsilon_\rho\rp$'
fields['eps_rho']['contour_name'] = r'$\varepsilon_\rho = 1$'
fields['eps_rho']['lims'] = np.array([-3, 1])
fields['eps_rho']['cmap'] = 'viridis'

fields['rp_over_ra'] = {}
fields['rp_over_ra']['array'] = np.log10(rp_over_ra)
fields['rp_over_ra']['tex_name'] = r'$\log_{10} \lp R_\text{SB}/R_a\rp$'
fields['rp_over_ra']['contour_name'] = r'$R_\text{SB} = R_a$'
fields['rp_over_ra']['lims'] = np.array([-2, 2])
fields['rp_over_ra']['cmap'] = 'PiYG'

for field in fields.keys():

  fig, ax = plt.subplots()
  
  ax.loglog()
  
  # Axes limits
  if plot_x == 'mass':
    ax.set_xlim(0.45, 0.8)
  else:
    ax.set_xlim(1e-1, np.amax(prof['r'] / cgs['RSUN']))
  ax.set_ylim(1e-2, 1e2)
  
  ax2 = ax.twiny()
  ax2.tick_params(which = 'both', direction = 'in')
  ax2.set_xlim(1e-1 * cgs['RSUN'] / cgs['AU'], np.amax(prof['r']) / cgs['AU'])
  ax2.set_xscale('log')
  ax2.set_xlabel(r'$a$ $\ls\unit{\astronomicalunit}\rs$', labelpad = 10)
  
  # Axes labels
  if plot_x == 'mass':
    ax.set_xlabel(r'$M_\text{enc}/M_\odot$')
  else:
    ax.set_xlabel(r'$a/R_\odot$')
  ax.set_ylabel(r'$M_\text{SB} / M_\text{Jup}$')
  
  # Axes ticks
  ax.minorticks_on()
  ax.tick_params(which = 'both', direction = 'in')
  
  # Plot
  x, y, z = rt.plot.pointsToPcolormeshArgs(plot_x_coord[min_idx:], mp / cgs['MJUP'], fields[field]['array'] )
  im = ax.pcolormesh(x, y, z, edgecolors = 'face', cmap = fields[field]['cmap'], vmin = fields[field]['lims'][0], vmax = fields[field]['lims'][1], rasterized = True)

  # eps_rho = 1 contour
  xC, yC = np.meshgrid(plot_x_coord, mp /cgs['MJUP'])
  cont4 = ax.contour(xC[:,:-500], yC[:,:-500], epsilon_rho[:-500,:].T, [1], colors = 'black', linestyles = 'dashed')
  ax.clabel(cont4, cont4.levels, fmt = r"$\varepsilon_\rho=1$", manual = [(.2, .2)], colors = 'black', fontsize = 'xx-large', use_clabeltext = True, inline_spacing = 8)
  
  # rp_over_ra = 1 contour
  cont = ax.contour(xC, yC, rp_over_ra.T, [1], colors = 'black', linestyles = 'dashed')
  ax.clabel(cont, cont.levels, fmt = r"$R_\text{SB} = R_a$", manual = [(3, 20)], colors = 'black', fontsize = 'xx-large', use_clabeltext = True, inline_spacing = 8)
  
  # a = Rt contour
  cont2 = ax.contour(xC, yC, a_over_rt.T, [1], colors = 'black', linestyles = 'dotted')
  ax.clabel(cont2, cont2.levels, fmt = r"$a = R_t$", manual = [(.5, 2)], colors = 'black', fontsize = 'xx-large', use_clabeltext = True, inline_spacing = 8)
  
  # rho v^2 = rho_p vesc^2
  cont3 = ax.contour(xC, yC, jiaspruit_f.T, [1], colors = 'black', linestyles = 'dotted')
  ax.clabel(cont3, cont3.levels, fmt = r"$f=1$", manual = [(.15, 2)], colors = 'black', fontsize = 'xx-large', use_clabeltext = True, inline_spacing = 8)

  if field == 'rp_over_ra':  
    ax.text(20, 20, 'Gravitational regime', fontsize = 'x-large', horizontalalignment = 'center')
    ax.text(20, 0.1, 'Geometrical regime', fontsize = 'x-large', horizontalalignment = 'center')
 
  # Colorbar
  cbar = plt.colorbar(im, pad = 0.01)
  cbar.set_label(fields[field]['tex_name'], rotation=90)
  cbar.ax.tick_params(axis='y', direction='in')
  
  plt.tight_layout()
  if field == 'eps_rho':
    side = 'right'
  else:
    side = 'left'
  plt.savefig('2_%s.pdf' % side, bbox_inches = 'tight')
  plt.close()
