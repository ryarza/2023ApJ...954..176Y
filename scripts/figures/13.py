import pathlib
import matplotlib
import matplotlib.figure
import matplotlib.style
import numpy as np
import rytools
import rytools.planets
import unyt

matplotlib.style.use(pathlib.Path.home() / '.config/matplotlib/style.mplstyle')


def au_of_rsun(x_in_rsun):
    return x_in_rsun * ((1 * unyt.r_sun).to(unyt.au) / unyt.au)


def rsun_of_au(x_in_au):
    return x_in_au * ((1 * unyt.au).to(unyt.r_sun) / unyt.r_sun)


def print_stats(planets: rytools.planets.ExoplanetArchiveData) -> None:
    print(f"Average planet stellar mass: {np.mean(planets.star_mass)}")
    frac = np.std(planets.star_mass) / np.mean(planets.star_mass)
    print(f"Fractional st dev planet stellar mass: {frac}")


def main() -> None:

    systems = rytools.planets.ExoplanetArchiveData(
        'data/planets.csv')
    bd_systems = rytools.planets.BrownDwarfData(
        'data/analytical/brown_dwarfs.csv')

    outcome_data = rytools.h5_to_dict('data/outcomes.h5')
    outcome_data['minimum_survive_conv_masses'] =\
        outcome_data['minimum_survive_conv_masses'] * unyt.m_jup
    outcome_data['minimum_ejection_masses'] =\
        outcome_data['minimum_ejection_masses'] * unyt.m_jup
    outcome_data['stellar_radii'] = outcome_data['stellar_radii'] * unyt.r_sun

    fig = matplotlib.figure.Figure()
    axis = fig.add_subplot()

    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylim(1e-1, 1e2)
    xlims = np.array([1e-2, 5]) * unyt.au
    axis.set_xlim(xlims)

    stellar_mass_mask = (systems.star_mass > 0.9 * unyt.m_sun) &\
        (systems.star_mass < 2 * unyt.m_sun)
    axis.scatter(
        systems.orbital_separation[stellar_mass_mask],
        systems.companion_mass[stellar_mass_mask],
        marker='x', color='black'
    )

    stellar_mass_mask = (bd_systems.star_mass > 0.9 * unyt.m_sun) &\
        (bd_systems.star_mass < 2 * unyt.m_sun)
    axis.scatter(
        bd_systems.orbital_separation[stellar_mass_mask],
        bd_systems.companion_mass[stellar_mass_mask],
        marker='x', color='black'
    )

    axis.plot(
        outcome_data['stellar_radii'].to(unyt.au),
        outcome_data['minimum_ejection_masses'].to(unyt.m_jup),
        color='black', ls='dashed'
    )

    axis.plot(
        outcome_data['stellar_radii'].to(unyt.au),
        outcome_data['minimum_survive_conv_masses'].to(unyt.m_jup),
        color='black', ls='dashed'
    )

    bbox = dict(facecolor='white', edgecolor='black', alpha=0.75)
    axis.text(
        0.3, 50, 'Envelope ejection',
        ha='center', fontsize='x-large', bbox=bbox
    )
    axis.text(
        .15, 7, "Destroyed below\nconvective region",
        ha='center', fontsize='x-large', bbox=bbox
    )
    axis.text(
        2.75e-2, 0.15, "Destroyed in\nconvective region",
        ha='center', fontsize='x-large', bbox=bbox
    )

    axis.axvline(
        outcome_data['stellar_radii'][-1].to(unyt.au),
        ls='dashed', color='grey'
    )
    axis.text(
        1.1 * outcome_data['stellar_radii'][-1].to(unyt.au),
        4 * axis.get_ylim()[0], 'Maximum RGB radius',
        rotation=-90, va='center', ha='center'
    )

    axis.set_xlabel(r'Orbital separation $\ls\unit{\astronomicalunit}\rs$')
    axis.set_ylabel(r'Companion mass $\ls M_\text{Jup}\rs$')

    ax2 = axis.twiny()
    ax2.set_xlim(xlims.to(unyt.r_sun))
    ax2.set_xscale('log')
    ax2.set_xlabel(r'Orbital separation $\ls R_\odot \rs$')

    fig.savefig('plots/outcomes.pdf')


main()
