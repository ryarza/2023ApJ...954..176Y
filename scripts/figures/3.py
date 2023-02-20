import argparse
import os
import pathlib
import matplotlib
import numpy as np
import rytools
import rytools.planets
import rytools.stars
import unyt

os.chdir(pathlib.Path(__file__).parent.parent)
matplotlib.style.use('~/.config/matplotlib/ricardo.mplstyle')

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help="Input profile", required=True)
args = parser.parse_args()

star = rytools.stars.MesaProfile(args.inputdata)

companion_masses = np.array([1, 10, 80]) * unyt.m_jup
companions = [rytools.planets.SubstellarBody(
    mass) for mass in companion_masses]

# Secondary axis with orbital separation


def m_of_a(separations):
    # if not isinstance(separations, collections.abc.Iterable):

    masses = np.empty_like(separations)
    for idx, separation in enumerate(separations):
        closest_idx, _ = rytools.nearest(
            star.r.to(unyt.r_sun).value,
            separation
        )
        masses[idx] = star.enclosed_mass[closest_idx].to(unyt.m_sun)
    # else
    #     out = prof.m_enc[rt.nearest(prof['r'], a * cgs['RSUN'])[0]]
    # return out
    return masses


def a_of_m(mval):
    # print(type(mval))
    # exit(0)
    # if type(mval) == np.ndarray:
    out = np.empty_like(mval)
    for out_idx in range(len(out)):
        closest_idx, _ = rytools.nearest(
            star.enclosed_mass.to(unyt.m_sun).value, mval[out_idx])
        out[out_idx] = star.r[closest_idx].to(unyt.r_sun)
    # else:
    #     out = prof['r'][rt.nearest(prof['menc'], mval * cgs['MSUN'])[0]]
    # return out / cgs['RSUN']
    return out


fig = matplotlib.figure.Figure()
axis = fig.add_subplot()

axis.plot(
    star.r.to(unyt.r_sun),
    - star.e_bind_outer_including_thermal.to(unyt.erg),
    # - star.e_bind_outer.to(unyt.erg),
    label=r'$E_\text{bind}^*$',
    color='black'
)

color_cycle = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']

for color, companion in zip(color_cycle, companions):
    ce = rytools.planets.AnalyticalCommonEnvelope(star, companion)

    axis.plot(
        ce.star.r[ce.companion_destruction_idx + 1:].to(unyt.r_sun),
        - ce.delta_e_orb[ce.companion_destruction_idx + 1:].to(unyt.erg),
        label=r'$\Delta E_\text{orb}^*\lp'
              f'{int(companion.mass.to(unyt.m_jup).value)}'
              r'M_\text{Jup}\rp$',
        color=color
    )

    axis.plot(
        ce.star.r[:ce.companion_destruction_idx + 1].to(unyt.r_sun),
        - ce.delta_e_orb[:ce.companion_destruction_idx + 1].to(unyt.erg),
        ls='dotted',
        color=color
    )

    # axis.plot(
    #     ce.star.r.to(unyt.r_sun),
    #     - ce.envelope_binding_energy.to(unyt.erg),
    #     # label=r'$\Delta E_\text{orb}\lp'
    #     #   f'{int(companion.mass.to(unyt.m_jup).value)}'
    #     #   r'M_\text{Jup}\rp$',
    #     color='red'
    # )

axis.axvline(star.core_radius().to(unyt.r_sun), ls='dashed', color='grey')
axis.text(0.0225, 1.5e45, 'Core-envelope boundary',
          rotation=90, va='center', ha='center')

axis.axhline(
    - star.e_bind_outer_including_thermal[star.core_idx()].to(unyt.erg),
    # - star.e_bind_outer[star.core_idx()].to(unyt.erg),
    ls='dashed',
    color='black'
)

axis.text(70, 3.5e46, 'Envelope binding energy',
          va='center', ha='center')

axis.set_xscale('log')
axis.set_yscale('log')
axis.set_xlim(1e-2, None)
axis.set_ylim(1e44, 1e49)
axis.set_xlabel(r'Orbital separation $\ls R_\odot \rs$')
axis.set_ylabel(r'$-\text{Energy}$ $\ls \unit{\erg} \rs$')

secax = axis.secondary_xaxis('top', functions=(m_of_a, a_of_m))

avals = np.array([1e-2, 1e-1, 1e0, 1e1, 1e2])
mvals = m_of_a(avals)
secax.set_xticks(mvals)
secax.set_xticklabels([f"{mass:.3f}" for mass in mvals])
secax.minorticks_off()
secax.set_xlabel(r'Enclosed mass $\ls M_\odot \rs$')

axis.legend(loc='upper right', ncol=2)

fig.savefig('plots/standard_energy_formalism.pdf')
