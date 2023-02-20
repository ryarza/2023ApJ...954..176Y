"""Plot analytical and numerical energy eposition rates."""
import os
import pathlib
import matplotlib
import matplotlib.figure
import matplotlib.style
import numpy as np
import rytools.planets
import unyt

os.chdir(pathlib.Path(__file__).parent.parent)

PATH = {}
PATH['GRID'] = pathlib.Path(
    'data/large_sims/inspirals/orbital_trajectories'
)


def check_consistency(
    sim: rytools.planets.EngulfmentIntegratorSimulation,
    companion_mass: unyt.array.unyt_quantity,
    use_drag_coefficients: int
) -> None:
    assert np.isclose(
        sim.companion.mass.to(unyt.m_jup).value,
        companion_mass.to(unyt.m_jup).value
    )
    assert use_drag_coefficients == sim.par['use_drag_coefficients']


def line_style(use_drag_coefficients: int) -> str:
    if use_drag_coefficients == 0:
        style = 'dashed'
    else:
        style = 'solid'
    return style


def main():
    """Do main calculation."""
    matplotlib.style.use('~/.config/matplotlib/style.mplstyle')
    for companion_mass in [1, 80] * unyt.mjup:
        for use_drag_coefficients in [0, 1]:
            # Load sim
            sim_path = PATH['GRID'] /\
                f"{use_drag_coefficients}" /\
                f"{int(companion_mass.to(unyt.mjup).value)}"
            print(f"Loading {sim_path}")
            sim = rytools.planets.EngulfmentIntegratorSimulation(sim_path)
            print(
                f"Inspiral time: {sim.time[:sim.destroyed_idx][-1].to(unyt.day)}")
            check_consistency(sim, companion_mass, use_drag_coefficients)
            # Plot
            fig = matplotlib.figure.Figure()
            axes = fig.add_subplot()
            if companion_mass == 1 * unyt.m_jup:
                axes.plot(
                    # self.x[:, 0].to(unyt.r_sun)[:self.destroyed_idx],
                    # self.x[:, 1].to(unyt.r_sun)[:self.destroyed_idx],
                    sim.position[:, 0].to(unyt.r_sun)[:sim.destroyed_idx],
                    sim.position[:, 1].to(unyt.r_sun)[:sim.destroyed_idx],
                    color='black',
                    lw=0.01,
                    # ls=line_style(use_drag_coefficients)
                )
                axes.plot(
                    sim.position[:, 0].to(unyt.r_sun)[sim.destroyed_idx:],
                    sim.position[:, 1].to(unyt.r_sun)[sim.destroyed_idx:],
                    color='black',
                    alpha=0.3,
                    lw=0.01,
                    # ls=line_style(use_drag_coefficients)
                )
            elif companion_mass == 80 * unyt.m_jup:
                max_idx = np.where(sim.r < 0.025 * sim.r[0])[0][0]
                axes.plot(
                    # self.x[:, 0].to(unyt.r_sun)[:self.destroyed_idx],
                    # self.x[:, 1].to(unyt.r_sun)[:self.destroyed_idx],
                    sim.position[:, 0].to(unyt.r_sun)[:max_idx],
                    sim.position[:, 1].to(unyt.r_sun)[:max_idx],
                    color='black',
                    lw=0.01,
                    # ls=line_style(use_drag_coefficients)
                )

            circle = matplotlib.patches.Circle(
                (0, 0),
                sim.prof.r[-1].to(unyt.r_sun),
                color='grey',
                fill=False,
                ls='dashed',
                label='Stellar surface'
            )
            axes.add_artist(circle)
            axes.set_xlabel(r'$x/R_\odot$')
            axes.set_ylabel(r'$y/R_\odot$')
            if companion_mass == 1 * unyt.m_jup:
                axes.set_xlim(-12, 12)
                axes.set_ylim(-12, 12)
                if not use_drag_coefficients:
                    axes.legend(loc=0)
            elif companion_mass == 80 * unyt.m_jup:
                axes.set_xlim(-110, 110)
                axes.set_ylim(-110, 110)
            axes.set_aspect(1)
            fig.savefig(
                "plots/inspiral"
                f"-{int(companion_mass.to(unyt.m_jup).value)}"
                f"-{use_drag_coefficients}.pdf")


main()
