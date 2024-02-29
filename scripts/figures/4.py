import argparse
import pathlib
import matplotlib
import matplotlib.figure
import matplotlib.style
import numpy as np
import rytools as rt
import scipy
import unyt


def soften_masses(minimum_mass: dict) -> tuple:
    x = minimum_mass['max_stellar_radii']
    y = minimum_mass['minimum_ejection_masses']
    y_f_1 = scipy.interpolate.interp1d(x, y, kind='quadratic')
    x_itp_1 = np.linspace(x[0], x[-1], 1000)
    y_itp_1 = y_f_1(x_itp_1)
    x_filt_1, y_filt_1 = scipy.signal.savgol_filter((x_itp_1, y_itp_1), 45, 3)
    print(np.amax(np.abs(x_filt_1 / x_itp_1 - 1)))
    return x_filt_1, y_filt_1


def au_of_rsun(x_in_rsun):
    return x_in_rsun * ((1 * unyt.r_sun).to(unyt.au) / unyt.au)


def rsun_of_au(x_in_au):
    return x_in_au * ((1 * unyt.au).to(unyt.r_sun) / unyt.r_sun)


parser = argparse.ArgumentParser()
parser.add_argument(
    "-i1", "--input1",
    help="Input file",
    required=False,
    default=pathlib.Path('data/rp_over_ra_at_destruction.h5')
)
parser.add_argument(
    "-i2", "--input2",
    help="Input file",
    required=False,
    default=pathlib.Path('data/minimum_ejection_mass.h5')
)
args = parser.parse_args()
matplotlib.style.use(pathlib.Path.home() / '.config/matplotlib/style.mplstyle')


def main() -> None:

    rp_over_ra = rt.h5_to_dict(args.input1)
    minimum_mass = rt.h5_to_dict(args.input2)

    assert (np.diff(rp_over_ra['stellar_radius']) > 0).all(),\
        np.where(np.diff(rp_over_ra['stellar_radius']) < 0)[0]

    rp_over_ra['stellar_radius'] = rp_over_ra['stellar_radius'] * unyt.r_sun
    rp_over_ra['planet_mass'] = rp_over_ra['planet_mass'] * unyt.m_jup

    idxs =\
        np.unravel_index(
            np.argmin(
                rp_over_ra['rp_over_ra'],
                axis=None
            ),
            rp_over_ra['rp_over_ra'].shape
        )

    print(rp_over_ra['stellar_radius'][idxs[0]])
    print(rp_over_ra['planet_mass'][idxs[1]])
    print(np.amin(rp_over_ra['rp_over_ra']))
    print(np.amax(rp_over_ra['rp_over_ra']))

    # x, y, z = rp_over_ra['stellar_radius'], rp_over_ra['planet_mass'], 2 * \
    # np.log10(2 * rp_over_ra['rp_over_ra'].T)
    x, y, z = rp_over_ra['stellar_radius'], rp_over_ra['planet_mass'], rp_over_ra['rp_over_ra'].T
    fig = matplotlib.figure.Figure()
    axis = fig.add_subplot()

    axis.set_xscale('log')
    # axis.set_yscale('log')

    # im = axis.pcolormesh(
    #     x, y, z,
    #     vmin=1, vmax=7,
    #     shading='nearest',
    #     rasterized=True
    # )

    im = axis.contourf(x, y, z, [0, 1, 2, 3, 4, 5, 6, 7], rasterized=True)

    cbar = fig.colorbar(im, pad=0.01)
    # cbar.set_label(r'$4\log_{10}\lp R_\text{SB} / R_a \rp$', rotation=90)
    cbar.set_label(
        r'$R_\text{SB} / R_a$ at the end of engulfment', rotation=90)

    axis.text(40, 55, 'Envelope ejection', color='white',
              ha='center', fontsize='xx-large')
    axis.text(15, 45, "Merger", color='white',
              ha='center', fontsize='xx-large')

    axis.set_xlim(left=10)
    # ax2 = ax.twiny()
    # secax = ax.secondary_xaxis('top', functions=(au_of_rsun, rsun_of_au))
    # secax.set_xlabel('angle [rad]')
    # ax2.tick_params(which = 'both', direction = 'in')
    # ax2.set_xlim(10 * cgs['RSUN'] / cgs['AU'], ax.get_xlim()[1] * cgs['RSUN'] / cgs['AU'])
    # ax2.set_xscale('log')
    # secax.set_xlabel(r'$a$ $\ls\unit{\astronomicalunit}\rs$', labelpad = 10)
    # secax.set_xscale('linear')

    smooth_radii, smooth_mass = soften_masses(minimum_mass)
    axis.plot(smooth_radii, smooth_mass, color='black')

    axis.set_ylim(10, 80)
    axis.set_xlabel(r'Stellar radius $\ls R_\odot\rs$')
    axis.set_ylabel(r'Companion mass $\ls M_\text{Jup}\rs$')

    ax2 = axis.twiny()
    ax2.set_xlim(axis.get_xlim())
    ax2.set_xscale('log')
    au_vals = np.array([0.05, 0.1, 0.2, 0.4, 0.6, 0.8])
    rsun_vals = rsun_of_au(au_vals)
    print(rsun_vals)
    ax2.set_xticks(rsun_vals)
    ax2.set_xticklabels(["%.2f" % au_vals[0]] + ["%.1f" %
                        i for i in au_vals[1:]])
    # ax2.minorticks_off()
    ax2.set_xlabel(
        r'Stellar radius $\ls\unit{\astronomicalunit}\rs$', labelpad=10)

    fig.savefig('plots/rp_over_ra_at_destruction.pdf')


main()
