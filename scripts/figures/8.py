"""Plot the parameter space explored by a engulfed substellar companion."""
import argparse
import pathlib
import matplotlib
import matplotlib.style
import numpy as np
import numpy.typing as npt
import rytools as rt
import rytools.units
import rytools.cewt
import unyt

matplotlib.style.use('~/.config/matplotlib/style.mplstyle')

GRID_PATH = pathlib.Path('data/large_sims/parameter_space/')
assert pathlib.Path(GRID_PATH).is_dir()


parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input",
    help="Input file",
    required=False,
    type=pathlib.Path,
    default=pathlib.Path('data/dimensionless_parameter_space.h5')
)

args = parser.parse_args()


def points_in_parameter_space(mass_ratio: npt.NDArray[np.float64],
                              eps_rho: npt.NDArray[np.float64],
                              rp_over_ra: npt.NDArray[np.float64]) -> tuple:
    # mask = np.where(
    #     (eps_rho < 1.5 * pow(q, 0.9)) &
    #     (eps_rho < 6 * pow(q, 0.425)) &
    #     (rp_over_ra < 0.5 * pow(q, - 2 / 3)) &
    #     (rp_over_ra < 0.1 * (1 + 1 / q)) &
    #     (rp_over_ra < 125 * pow(eps_rho, 0.9))
    # )
    mask = (rp_over_ra < 0.5 * pow(mass_ratio, - 2 / 3)) &\
        (rp_over_ra < 0.1 * (1 + 1 / mass_ratio)) &\
        (eps_rho < 6 * pow(mass_ratio, 0.425)) &\
        (eps_rho > 1.5 * pow(mass_ratio, 0.9)) &\
        (rp_over_ra < 125 * pow(eps_rho, 0.9)) &\
        (rp_over_ra > 4.8e-3 * pow(eps_rho, -0.9)) &\
        (rp_over_ra < 4.8 * pow(eps_rho, -1.5)
         )
    return mass_ratio[mask], eps_rho[mask], rp_over_ra[mask]


def main():
    """Do main calculation."""
    print(f"Plotting {args.input}")

    d = rt.h5_to_dict(args.input)
    d['sb_masses'] = d['sb_masses'] * unyt.m_jup
    d['max_stellar_radii'] = d['max_stellar_radii'] * unyt.r_sun
    print("Loading simulation grid...")
    grid = rytools.cewt.SimulationGrid(GRID_PATH)
    print("Loaded simulation grid...")
    q = grid["sim_q"]
    eps_rho = grid.eps_rho
    rp_over_ra = grid["sim_radius"]

    prop_cycle = matplotlib.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    linestyles = ['solid', 'dotted', 'dashed', 'dashdot']

    q, eps_rho, rp_over_ra =\
        points_in_parameter_space(q, eps_rho, rp_over_ra)

    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot()
    for i in range(d["q"].shape[0]):
        for j in range(d["q"].shape[1]):
            ax.plot(d["q"][i, j], d["eps_rho"][i, j],
                    lw=0.75, ls=linestyles[i], color=colors[j])
            ax.plot(d["q"][i, j][0], d["eps_rho"][i, j][0],
                    marker='^', color=colors[j], markersize=2.5)

    ax.scatter(q, eps_rho, marker="s", s=12,
               facecolor='none', edgecolor='black', lw=0.75)

    ax.loglog()
#    ax.legend(loc=0)
    ax.set_xlabel("$q$")
    ax.autoscale(enable=False)
    ax.fill_between([1 / 9, ax.get_xlim()[1]], *
                    ax.get_ylim(), color='gray', alpha=0.3, lw=0)
    ax.set_ylabel(r"$\varepsilon_\rho$")
    fig.savefig("plots/q_vs_eps_rho.pdf")

    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot()
    ax.loglog()
    custom_lines = [matplotlib.lines.Line2D(
        [0], [0], color='black', ls=ls) for ls in linestyles]
    for i in range(d["q"].shape[0]):
        for j in range(d["q"].shape[1]):
            ax.plot(d["q"][i, j], d["rp_over_ra"]
                    [i, j], color=colors[j], lw=0.75, ls=linestyles[i])
            ax.plot(d["q"][i, j][0], d["rp_over_ra"][i, j][0],
                    marker='^', color=colors[j], markersize=2.5)

    ax.scatter(q, rp_over_ra, marker="s", s=12,
               facecolor='none', edgecolor='black', lw=0.75)
    ax.legend(custom_lines, [
              rf'$R_\star={int(r)}R_\odot$' for r in d["max_stellar_radii"]])
    qs = np.logspace(-3, np.log10(1 / 9))
#    y = 0.5 * pow(qs, -2 / 3)
#    ax.plot(qs, y, ls="dashed")
#    ax.plot(qs, 0.1 * (1 + 1 / qs), ls="dashed", color='black')
#    ax.vlines(1 / 9, ymin=1e-3, ymax=1, ls='dashed', color='black')
    ax.set_ylim(1e-3, 1e2)
    ax.set_xlabel("$q$")
    ax.set_ylabel(r"$R_\text{SB}/R_a$")
    ax.autoscale(enable=False)
    ax.text(1e-1, 30, 'Region not suitable for\nwind tunnel simulations',
            ha='center', va='center')
    ax.fill_between(qs, 0.1 * (1 + 1 / qs), ax.get_ylim()
                    [1], color='gray', alpha=0.3, lw=0)
    ax.fill_between([1 / 9, ax.get_xlim()[1]], *
                    ax.get_ylim(), color='gray', alpha=0.3, lw=0)
    ax.fill_between(
        [*ax.get_xlim()],
        [ax.get_ylim()[0]] * 2,
        [0.2] * 2,
        color='blue',
        alpha=0.1,
        lw=0
    )
#    ax.text(6e-3, 6e-3, "Gravitational regime.\n"
#            r"Drag independent of $R_\text{SB}/R_a$",
#            ha='center', fontsize='x-large')
    fig.savefig("plots/q_vs_rp_over_ra.pdf")

    fig = matplotlib.figure.Figure()
    ax = fig.add_subplot()
    for j in range(d["q"].shape[1]):
        for i in range(d["q"].shape[0]):
            ax.plot(d["eps_rho"][i, j], d["rp_over_ra"]
                    [i, j], color=colors[j], lw=0.75, ls=linestyles[i])
            ax.plot(d["eps_rho"][i, j][0], d["rp_over_ra"][i, j][0],
                    marker='^', color=colors[j], markersize=2.5)
        ax.text(d["eps_rho"][-1, j][0] * 1.4, d["rp_over_ra"][-1, j][0],
                f"${round(float(d['sb_masses'][j].value))}"
                r"M_\text{Jup}$", color=colors[j], ha='center', va='center')

    ax.scatter(eps_rho, rp_over_ra, marker="s", s=12,
               facecolor='none', edgecolor='black', lw=0.75)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.autoscale(enable=False)

    ax.fill_between(
        [*ax.get_xlim()],
        [ax.get_ylim()[0]] * 2,
        [0.2] * 2,
        color='blue',
        alpha=0.1,
        lw=0
    )

    ax.text(
        2e-2,
        1e-2,
        'Gravitational regime\n'
        r'Weak dependence on $R_\text{SB}/R_a$',
        ha='center',
        va='center'
    )
#    ax.fill_between(np.linspace(*ax.get_xlim()),
#                    ax.get_ylim()[0], 0.3, color='blue', alpha=0.3)
#    ax.fill_between(np.linspace(*ax.get_xlim()), 5,
#                    ax.get_ylim()[1], color='green', alpha=0.3)

    ax.set_xlabel(r"$\varepsilon_\rho$")
    ax.set_ylabel(r"$R_\text{SB}/R_a$")
    fig.savefig("plots/eps_rho_vs_rp_over_ra.pdf")


main()
