"""Gravitational and geometrical regimes in a simulation grid."""
import itertools
import os
import pathlib
import matplotlib
import matplotlib.style
import matplotlib.figure
import numpy as np
import rytools.cewt

matplotlib.style.use('~/.config/matplotlib/style.mplstyle')
os.chdir(pathlib.Path(__file__).parent.parent)


def main() -> None:
    """Do main calculation."""
    grid_path = pathlib.Path('data/large_sims/parameter_space/q-2.15443e-02')
    # grid_path = pathlib.Path('large_sims/parameter_space/q-4.64159e-03')
    # grid_path = pathlib.Path('large_sims/parameter_space/q-1.00000e-01')
    grid = rytools.cewt.SimulationGrid(grid_path)
    markers = itertools.cycle(('o', 'x', 's', '^'))
    # grid = pyflash.SimulationGrid('large_sims/parameter_space/')

    # Choose only simulations that have ram pressure computed
    paths = [sim.path for sim in grid if
             (sim.path / "ram_pressure.h5").is_file() and sim.complete]

    grid = rytools.cewt.SimulationGrid(paths)

    # for sim in grid.sims:
    #     print(sim.path)
    #     cewt.load_cg(sim.path, 1.6 * max(1, sim.par["sim_radius"]))

    ram_coeff = grid.c_ram_av(15)
    grav_coeff = grid.c_grav_av(15, 1.6)
    rp_over_ras = grid['sim_radius']
    eps_rhos = grid.eps_rho

    my_eps_rhos = np.logspace(-1, 0, 4)

    fig = matplotlib.figure.Figure()
    axes = fig.add_subplot()
    axes.set_xscale('log')
    for eps in my_eps_rhos:
        mask = np.where(np.isclose(eps_rhos, eps))
        idxs = np.argsort(rp_over_ras[mask])
        axes.plot(
            rp_over_ras[mask][idxs],
            ram_coeff[mask][idxs] *
            np.minimum(1, pow(grid['sim_radius'][mask][idxs], 2)),
            label=fr"$\varepsilon_\rho={eps:.2f}$",
            marker=next(markers)
        )
    axes.text(0.4, -0.1, 'Gravitational regime\nNegligible ram pressure drag',
              ha='center', va='center', fontsize='x-large')
    axes.text(3, 0.25, 'Ram pressure regime\n'
              r'Approaching subsonic flow',
              ha='center', va='center', fontsize='x-large')
    axes.set_xlabel(r'$R_\text{SB}/R_a$')
    axes.set_ylabel(r'$C_p\min\lp1,\ls R_\text{SB}/R_a\rs^{2}\rp$')
    axes.set_ylim(-0.2, 1)
    axes.set_xlim(1e-1, 1e1)
    axes.legend(loc=0)
    fig.savefig("plots/cp.pdf")

    fig = matplotlib.figure.Figure()
    axes = fig.add_subplot()
    for eps in my_eps_rhos:
        mask = np.where(np.isclose(eps_rhos, eps))
        idxs = np.argsort(rp_over_ras[mask])
        axes.plot(rp_over_ras[mask][idxs],
                  grav_coeff[mask][idxs] *
                  np.minimum(1, pow(rp_over_ras[mask][idxs], -2)),
                  label=fr"$\varepsilon_\rho={eps:.2f}$",
                  marker=next(markers))
    axes.set_xlabel(r'$R_\text{SB}/R_a$')
    axes.set_ylabel(r'$C_g\min\lp1,\ls R_\text{SB}/R_a\rs^{-2}\rp$')
    axes.text(3, 0.15, 'Ram pressure regime\nNegligible gravitational drag',
              ha='center', va='center', fontsize='x-large')
    axes.text(0.4, 1.35, 'Gravitational regime\n'
              r'Weak dependence on $R_\text{SB}/R_a$',
              ha='center', va='center', fontsize='x-large')
    axes.set_xscale('log')
    axes.set_xlim(1e-1, 1e1)
    axes.set_ylim(-0.5, 1.5)
    fig.savefig("plots/cg.pdf")


main()
