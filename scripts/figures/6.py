"""Plot a grid of density slices at varying density gradients and sizes."""
import pathlib
import typing
import numpy as np
import matplotlib.style
import yt
import rytools as rt

matplotlib.style.use('~/.config/matplotlib/ricardo.mplstyle')


def plot_variables() -> typing.Tuple:
    """Values for the independent variables at which to plot slices."""
    eps_rhos = np.logspace(-1, 0, 4)[[0, 2, -1]]
    rp_over_ras = np.logspace(np.log10(0.3), 0, 7)
    rp_over_ras = pow(10, rt.expand_linspace(
        np.log10(rp_over_ras), 1, 6))
    rp_over_ras = rp_over_ras[[0, 6, -1]]
    return eps_rhos, rp_over_ras


def yt_dataset_at_time(sim_path: pathlib.Path, plot_time: float):
    """Return the yt dataset with time closest to provided time."""
    plot_file = 100
    checkpoint_path = sim_path / f'acc_hdf5_chk_{str(plot_file).zfill(4)}'
    dataset = yt.load(checkpoint_path)
    if not np.isclose(float(dataset.current_time), plot_time, rtol=1.e-1):
        plot_file = round(plot_file * plot_time / float(dataset.current_time))
        checkpoint_path = sim_path / f'acc_hdf5_chk_{str(plot_file).zfill(4)}'
        dataset = yt.load(checkpoint_path)
        assert np.isclose(float(dataset.current_time), plot_time, rtol=1.e-1),\
            f"Couldn't find checkpoint file at requested time {plot_time}. "\
            f"Closest was {dataset.current_time}"

    print(f"  Plotting {checkpoint_path} with time {dataset.current_time}")
    return dataset


def plot_panel(axis,
               edge: typing.List[bool],
               rp_over_ra: float,
               eps_rho: float,
               yt_dataset):
    """Plot a single panel in the grid."""
    cax = axis
    cax.set_aspect(1)
    # If plot belongs to bottom row
    if edge[0]:
        cax.set_xlabel(r'$x/\max\lp R_\text{SB}, R_a\rp$')
    else:
        cax.xaxis.set_ticklabels([])
    # If panel belongs to leftmost column
    if edge[1]:
        cax.set_ylabel(r'$y/\max\lp R_\text{SB}, R_a\rp$')
    else:
        cax.yaxis.set_ticklabels([])

    l_domain_half = 4 * max(1, rp_over_ra)
    x_0, y_0 = -l_domain_half, -l_domain_half
    x_1, y_1 = l_domain_half, l_domain_half
    pixels_x = 1024
    pixels_y = int(pixels_x * (y_1 - y_0) / (x_1 - x_0))
    x = np.linspace(x_0, x_1, pixels_x)
    y = np.linspace(y_0, y_1, pixels_y)

    slc = yt_dataset.slice("z", 0)

    frb = yt.visualization.fixed_resolution.FixedResolutionBuffer(
        slc, (x_0, x_1, y_0, y_1), (pixels_x, pixels_y))

    rho = np.array(frb["gas", "density"]).T
    bdry_var = np.array(frb['flash', 'bdry']).T
    vel_x = np.array(frb['gas', 'velocity_x']).T
    vel_y = np.array(frb['gas', 'velocity_y']).T
    vel_z = np.array(frb['gas', 'velocity_z']).T

    rho[bdry_var > 0] = np.nan

#               print(e)
#               rho = np.ones((len(x), len(y)))
#               vel_x = 1 * rho
#               vel_y = 0 * rho
#               vel_z = 0 * rho
#               bdry_var = -1 * rho

    vel_abs = np.sqrt(vel_x * vel_x + vel_y * vel_y + vel_z * vel_z)
    vel_abs[bdry_var > 0] = 0

    x_scaled = x / max(1, rp_over_ra)
    y_scaled = y / max(1, rp_over_ra)
    x_pcm, y_pcm = np.meshgrid(x_scaled, y_scaled)
    # Background color
    cax.set_facecolor("black")

    cax.set_xlim(*x_scaled[[0, -1]])
    cax.set_ylim(*y_scaled[[0, -1]])
    pcm = cax.pcolormesh(x_pcm, y_pcm, np.log10(rho).T, shading='nearest',
                         vmin=-2, vmax=1, rasterized=True, cmap='inferno')
    lw = 2 * vel_abs / np.amax(vel_abs)
    cax.streamplot(x_scaled, y_scaled, vel_x.T, vel_y.T,
                   linewidth=lw.T, color='black')

    cax.plot([0], [0], visible=False,
             label=fr'$\varepsilon_\rho ={eps_rho:.2f}$,'
             r' $R_\text{SB}/R_a='f'{rp_over_ra:.2f}$')
    leg = cax.legend(handlelength=0, handletextpad=0, loc='upper left')
    for item in leg.legendHandles:
        item.set_visible(False)
    return cax, pcm


def main() -> None:
    """Do main calculation."""
    rt.set_work_dir(pathlib.Path(__file__))
    grid_path = pathlib.Path('large_sims/parameter_space/q-2.15443e-02')
    # grid_path = pathlib.Path('large_sims/parameter_space/q-4.64159e-03')
    assert grid_path.is_dir()

    eps_rhos, rp_over_ras = plot_variables()

    # fig, ax = plt.subplots(nrows = len(rp_over_ras), ncols = len(eps_rhos),
    # constrained_layout = True, figsize = (5.45 * 3, 4.8 * 3))
    # fig = matplotlib.figure.Figure(constrained_layout=True,
    #                                figsize=(4 * 3, 3.5 * 3)
    #                                )
    fig = matplotlib.figure.Figure(constrained_layout=True,
                                   figsize=(4 * 3 * 1.01, 3.5 * 3 * 1.01)
                                   )

    axes = fig.subplots(ncols=len(rp_over_ras), nrows=len(eps_rhos))

    for idx_rp_over_ra, rp_over_ra in enumerate(rp_over_ras):
        print(f"Plotting Rp / Ra = {float(rp_over_ra):.2f}")
        for idx_eps_rho, eps_rho in enumerate(eps_rhos):
            print(f"  Plotting eps_rho = {float(eps_rho):.2f}")
            sim_path = grid_path / \
                f"rp_over_ra-{rp_over_ra:.5e}" / f"eps_rho-{eps_rho:.5e}"
            yt_dataset = yt_dataset_at_time(sim_path, 20 * max(1, rp_over_ra))
            axis = axes[idx_eps_rho, idx_rp_over_ra]
            axes[idx_eps_rho, idx_rp_over_ra], pcm =\
                plot_panel(axis,
                           [
                               idx_eps_rho == len(eps_rhos) - 1,
                               idx_rp_over_ra == 0
                           ],
                           rp_over_ra, eps_rho, yt_dataset)

    cbar = fig.colorbar(pcm, ax=axes.ravel().tolist(), aspect=30,
                        label=r'$\log_{10} \lp \rho / \rho_\infty \rp$',
                        pad=0.01)
    cbar.ax.tick_params(which='both', direction='in')

    fig.savefig('plots/grid_1.pdf')


main()
