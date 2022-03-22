import numpy as np
import matplotlib.pyplot as plt
import yt
import rytools as rt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Directory with simulations", required = True)
args = parser.parse_args()

rt.plot.plottex(fontsize = 'xx-large')

eps_rhos = ['0.1', '0.5', '1']
rsb_over_ras = ['0.3', '0.7', '1']
plotTime = 25
plotFile = 100

fig, ax = plt.subplots(ncols = len(rsb_over_ras), nrows = len(eps_rhos), constrained_layout = True, figsize = (4 * 3, 3.5 * 3))

for idx_rsb_over_ra, rsb_over_ra in enumerate(rsb_over_ras):
  print("Plotting simulations with Rsb / Ra = %.2f" % float(rsb_over_ra))
  for idx_eps_rho, eps_rho in enumerate(eps_rhos):
    print("  Plotting eps_rho = %.2f" % float(eps_rho))
    cax = ax[idx_eps_rho, idx_rsb_over_ra]
    cax.set_aspect(1)
    cax.tick_params(which = 'both', direction = 'in')
    if idx_eps_rho < len(eps_rhos) - 1: cax.xaxis.set_ticklabels([])
    else: cax.set_xlabel(r'$x/R_a$')
    if idx_rsb_over_ra > 0: cax.yaxis.set_ticklabels([])
    else: cax.set_ylabel(r'$y/R_a$')
    x_0, y_0, x_1, y_1 = -4, -4, 4, 4
    pixels_x = 256
    pixels_y = int(pixels_x * ( y_1 - y_0 ) / ( x_1 - x_0 ) )
    x = np.linspace(x_0, x_1, pixels_x)
    y = np.linspace(y_0, y_1, pixels_y)

    try:
      ds = yt.load(args.inputdata + eps_rho + '/' + rsb_over_ra + '/acc_hdf5_chk_' + str(plotFile).zfill(4))
      if np.abs(float(ds.current_time) / plotTime - 1) > 1e-1:
        plotFile = round(plotTime * plotTime / float(ds.current_time))
        ds = yt.load(args.inputdata + eps_rho + '/' + rsb_over_ra + '/acc_hdf5_chk_' + str(plotFile).zfill(4))
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

    cax.plot([0], [0], visible = False, label = r'$\varepsilon_\rho = ' + eps_rho + r', R_\text{SB}/R_a=' + rsb_over_ra + r'$')
    leg = cax.legend(handlelength=0, handletextpad=0, loc = 'upper left')
    for item in leg.legendHandles:
      item.set_visible(False)

cbar = fig.colorbar(pcm, ax=ax.ravel().tolist(), aspect = 30, label = r'$\log_{10} \lp \rho / \rho_\infty \rp$', shrink = 0.99)
cbar.ax.tick_params(which = 'both', direction = 'in')

print("Saving figure...")
plt.savefig('6.pdf', bbox_inches = 'tight')
plt.close()
