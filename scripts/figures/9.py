import numpy as np
import matplotlib.pyplot as plt
import yt
import rytools as rt
import argparse

#Input argument: directory to analyze
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Directory with simulations", required = True)
args = parser.parse_args()

rt.plot.plottex(fontsize = 17.28 * 0.75 * (16.35 * 13.7) / 12 / 10.5 )

qs = np.logspace(-3, -1, 4)
rp_over_ras = np.logspace(np.log10(0.3), 0, 4)
eps_rhos = np.logspace(-1, 0, 4)

rp_over_ra = rp_over_ras[0]

plotTime = 25
plotFile = 100
resolution = 1024

fig, ax = plt.subplots(ncols = len(qs), nrows = len(eps_rhos), constrained_layout = True, figsize = (16.35, 13.7))

for idx_q, q in enumerate(qs):
  print("Plotting simulations with q = %.2e" % float(q))
  for idx_eps_rho, eps_rho in enumerate(eps_rhos):
    print("  Plotting eps_rho = %.2f" % float(eps_rho))
    simDirectory = args.inputdata + 'q-%.5e/eps_rho-%.5e/rsb_over_ra-%.5e/' % ( q, eps_rho, rp_over_ra )
    cax = ax[idx_eps_rho, idx_q]
    cax.set_aspect(1)
    cax.tick_params(which = 'both', direction = 'in')
    if idx_eps_rho < len(eps_rhos) - 1: cax.xaxis.set_ticklabels([])
    else: cax.set_xlabel(r'$x/R_a$')
    if idx_q > 0: cax.yaxis.set_ticklabels([])
    else: cax.set_ylabel(r'$y/R_a$')
    x_0, y_0, x_1, y_1 = -4, -4, 4, 4
    pixels_x = resolution
    pixels_y = int(pixels_x * ( y_1 - y_0 ) / ( x_1 - x_0 ) )
    x = np.linspace(x_0, x_1, pixels_x)
    y = np.linspace(y_0, y_1, pixels_y)

    try:
      ds = yt.load(simDirectory + 'acc_hdf5_chk_' + str(plotFile).zfill(4))
      if np.abs(float(ds.current_time) / plotTime - 1) > 1e-1:
        plotFile = round(plotTime * plotTime / float(ds.current_time))
        ds = yt.load(simDirectory + 'acc_hdf5_chk_' + str(plotFile).zfill(4))
        assert np.abs(float(ds.current_time) / plotTime - 1) < 1e-2, "Couldn't find checkpoint file at requested time"
      
      slc = ds.slice("z", 0)

      frb = yt.visualization.fixed_resolution.FixedResolutionBuffer(slc, (x_0, x_1, y_0, y_1), (pixels_x, pixels_y))

      rho = np.array(frb["gas", "density"]).T
      bdry_var = np.array(frb['flash', 'bdry']).T
      vel_x = np.array(frb['gas', 'velocity_x']).T
      vel_y = np.array(frb['gas', 'velocity_y']).T
      vel_z = np.array(frb['gas', 'velocity_z']).T

      rho[bdry_var > 0] = np.nan

    except Exception as e:
      print(e)
      rho = np.ones((len(x), len(y)))
      vel_x = 1 * rho
      vel_y = 0 * rho
      vel_z = 0 * rho
      bdry_var = -1 * rho

    vel_abs = np.sqrt(vel_x * vel_x + vel_y * vel_y + vel_z * vel_z)
    vel_abs[bdry_var > 0] = 0

    x_pcm, y_pcm, rho_pcm = rt.plot.pointsToPcolormeshArgs(x, y, rho)
    # Background color
    cax.set_facecolor("black")

    cax.set_xlim(x_0, x_1)
    cax.set_ylim(y_0, y_1)
    pcm = cax.pcolormesh(x_pcm, y_pcm, np.log10(rho_pcm), vmin = -2, vmax = 1, rasterized = True, cmap = 'inferno')
    lw = 2 * vel_abs / np.amax(vel_abs)
    cax.streamplot(x, y, vel_x.T, vel_y.T, linewidth = lw.T, color = 'black')

    q_n2 = np.floor(np.log10(q))
    q_n1 = q / pow(10, q_n2)
    cax.plot([0], [0], visible = False, label = r'$\varepsilon_\rho = %.2f, q = %.1f\times10^{%i}$' % ( eps_rho, q_n1, q_n2) )
    leg = cax.legend(handlelength=0, handletextpad=0, loc = 'upper left', fontsize = 'xx-large')
    for item in leg.legendHandles:
      item.set_visible(False)

cbar = fig.colorbar(pcm, ax=ax.ravel().tolist(), aspect = 30, label = r'$\log_{10} \lp \rho / \rho_\infty \rp$', shrink = 0.99)
cbar.ax.tick_params(which = 'both', direction = 'in')

plt.savefig('8.pdf', bbox_inches = 'tight')
plt.close()
