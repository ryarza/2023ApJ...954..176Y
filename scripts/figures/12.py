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
    'data/large_sims/inspirals/compare_energy_deposition_rates/0')


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


def load_luminosities(
    sim: rytools.planets.EngulfmentIntegratorSimulation
) -> tuple[list, list]:
    times = np.loadtxt(
        sim.path / "luminosity_time.txt"
    ) * unyt.s
    luminosities = np.loadtxt(
        sim.path / "luminosity_luminosity.txt"
    ) * unyt.erg / unyt.s
    return times, luminosities


def compute_luminosities(
    sim: rytools.planets.EngulfmentIntegratorSimulation
) -> None:
    time, luminosity = sim.luminosity
    np.savetxt(
        sim.path / "luminosity_time.txt",
        time.to(unyt.s).value
    )
    np.savetxt(
        sim.path / "luminosity_luminosity.txt",
        luminosity.to(unyt.erg / unyt.s).value
    )


def main():
    """Do main calculation."""
    mosaic = [['A', 'C'], ['B', 'C']]
    matplotlib.style.use('~/.config/matplotlib/style.mplstyle')
    use_cached = False
    # fig = matplotlib.figure.Figure(figsize=(10, 4.8))
    fig = matplotlib.figure.Figure(figsize=(10 * 1.25, 4.8 * 1.25))
    axes = fig.subplot_mosaic(mosaic=mosaic)
    colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']
    # axes = fig.add_subplot()
    for companion_idx, companion_mass in enumerate([1, 4, 20, 80] * unyt.mjup):
        for use_drag_coefficients in [1]:
            sim_path = PATH['GRID'] /\
                f"{use_drag_coefficients}" /\
                f"{int(companion_mass.to(unyt.mjup).value)}"
            print(f"Loading {sim_path}")
            sim = rytools.planets.EngulfmentIntegratorSimulation(sim_path)
            print(f"Destroyed at {sim.r[sim.destroyed_idx].to(unyt.r_sun)}")
            print(f"{len(sim.r)} times")
            print(
                f"tkh at destruction: {sim.t_kh[sim.destroyed_idx].to(unyt.year)}")
            check_consistency(sim, companion_mass, use_drag_coefficients)
            if not use_cached:
                compute_luminosities(sim)
            time, luminosity = load_luminosities(sim)
            axes['C'].plot(
                time.to(unyt.year),
                1 + luminosity / sim.prof.luminosity[-1],
                label=f'${int(companion_mass.to(unyt.m_jup).value)}'
                r'M_\text{Jup}$'
            )
            # x, y = scipy.signal.savgol_filter(
            #     (
            #         sim.r.to(unyt.r_sun),
            #         np.abs(sim.orbital_decay_timescale.to(unyt.yr))
            #     ),
            #     500,
            #     polyorder=2
            # )
            axes['A'].plot(
                sim.r.to(unyt.r_sun)[::10],
                np.abs(sim.orbital_decay_timescale.to(unyt.day))[::10]
            )
            axes['B'].plot(
                sim.r.to(unyt.r_sun)[::10],
                sim.energy_deposition_rate.to(unyt.erg / unyt.s)[::10],
                ls='dashed',
                color=colors[companion_idx]
            )
            axes['B'].plot(
                sim.r.to(unyt.r_sun)[::10],
                luminosity.to(unyt.erg / unyt.s)[:len(sim.time):10],
                color=colors[companion_idx]
            )
            # x = sim.r.to(unyt.r_sun)
            # y = np.convolve(np.abs(sim.orbital_decay_timescale.to(
            # unyt.yr)), np.ones(100)/100, mode='valid')
            # y = scipy.ndimage.uniform_filter1d(
            #     np.abs(sim.orbital_decay_timescale.to(unyt.yr)),
            #     size=2000)
            # axes['A'].plot(
            #     x,
            #     y,
            #     lw=0.01
            # )

    axes['B'].sharex(axes['A'])
    axes['A'].plot(
        sim.prof.r.to(unyt.r_sun),
        sim.prof.t_kh.to(unyt.day),
        color='black',
        ls='dashed',
        label='Energy transport'
    )
    axes['A'].plot([0], [0], label=r'Orbital decay', ls='solid', color='black')
    axes['A'].set_ylabel(r'Timescale $\ls\unit{\day}\rs$')
    axes['A'].legend(loc=0)
    axes['A'].set_yscale('log')

    axes['A'].tick_params('x', labelbottom=False)
    axes['A'].set_ylim(1e-3, 1e4)
    axes['A'].set_xlim(0, 10)

    axes['B'].set_yscale('log')
    axes['B'].plot([0], [0], label=r'Additional luminosity',
                   ls='solid', color='black')
    axes['B'].plot([0], [0], label=r'Energy deposition rate',
                   ls='dashed', color='black')
    axes['B'].legend(loc=0)
    axes['B'].set_ylabel(r'$\ls \unit{\erg\per\second}\rs$')
    axes['B'].set_xlabel(r'Radial coordinate $\ls R_\odot\rs$')

    # axes.loglog()
    axes['C'].set_xscale('log')
    axes['C'].set_yscale('log')
    axes['C'].legend(loc=0)
    axes['C'].set_ylim(1, 1e5)
    axes['C'].set_xlim(1e-3, 1e4)
    axes['C'].set_xlabel(r'Time $\ls\unit{\year}\rs$')
    axes['C'].set_ylabel(r'Luminosity $\ls\text{stellar luminosity}\rs$')
    fig.savefig('plots/luminosity.pdf')


main()
