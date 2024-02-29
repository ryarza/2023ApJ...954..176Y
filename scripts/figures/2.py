import argparse
import pathlib
import warnings
import mesa_reader
import numpy as np
import rytools
import rytools.planets
import rytools.stars
import matplotlib
import matplotlib.figure
import matplotlib.style
import unyt

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", "--input",
    help="Input folder with MESA files",
    type=pathlib.Path,
    default=pathlib.Path(
        '/data/groups/ramirez-ruiz/rcastroy/mesa/'
        'grid_sun_like/1.00000e+00/LOGS_to_start_he_core_flash/'
    )
)
parser.add_argument(
    "-r", "--radius",
    help="Desired stellar radius (in rsun)",
    type=int,
    default=10
)

args = parser.parse_args()

matplotlib.style.use(pathlib.Path.home() / '.config/matplotlib/style.mplstyle')


def get_star() -> rytools.stars.MesaProfile:
    logdir = mesa_reader.MesaLogDir(args.input)
    profile_number = rytools.nearest(logdir.history.R, args.radius)[0] + 1
    profile_path = args.input / f"profile{profile_number}.data"
    star = rytools.stars.MesaProfile(profile_path)
    if not np.isclose(star.r[-1].to(unyt.r_sun).value, args.radius, rtol=1e-2):
        warnings.warn('Close stellar radius not found')
    return star


def main() -> None:
    star = get_star()

    print(f"Radius: {star.r[-1]}")
    print(f"Luminosity: {star.luminosity[-1]}")

    n_points = 200
    min_idx = 0
    # Range of planet masses (in Jupiter masses)
    companion_masses = np.logspace(-1, 2, n_points) * unyt.m_jup
    # Corresponding range of radii (in Jupiter radii)
    companions = [rytools.planets.SubstellarBody(mass)
                  for mass in companion_masses]

    ram_disruption_separations = np.empty(len(companions))
    tidal_disruption_separations = np.empty(len(companions))
    destruction_separations = np.empty(len(companions))
    epsilon_rho = np.empty((star.n_cells, len(companions)))
    rp_over_ra = np.empty((star.n_cells, len(companions)))
    for idx, companion in enumerate(companions):
        ce = rytools.planets.AnalyticalCommonEnvelope(star, companion)
        ram_disruption_separations[idx] = star.r[ce.destruction_idx_ram_pressure]
        tidal_disruption_separations[idx] = star.r[ce.destruction_idx_tides]
        destruction_separations[idx] = star.r[ce.companion_destruction_idx]
        epsilon_rho[:, idx] = ce.epsilon_rho
        rp_over_ra[:, idx] = ce.rsb_over_ra

    fields = {}
    fields['eps_rho'] = {}
    fields['eps_rho']['array'] = np.log10(epsilon_rho)
    fields['eps_rho']['tex_name'] = r'$\log_{10}\lp\varepsilon_\rho\rp$'
    fields['eps_rho']['contour_name'] = r'$\varepsilon_\rho = 1$'
    fields['eps_rho']['lims'] = np.array([-3, 1])
    fields['eps_rho']['cmap'] = 'viridis'

    fields['rp_over_ra'] = {}
    fields['rp_over_ra']['array'] = np.log10(rp_over_ra)
    fields['rp_over_ra']['tex_name'] = r'$\log_{10} \lp R_\text{SB}/R_a\rp$'
    fields['rp_over_ra']['contour_name'] = r'$R_\text{SB} = R_a$'
    fields['rp_over_ra']['lims'] = np.array([-2, 2])
    fields['rp_over_ra']['cmap'] = 'PiYG'

    # fields['mach'] = {}
    # fields['mach']['array'] = mach
    # fields['mach']['tex_name'] = r'Mach number'
    # fields['mach']['contour_name'] = r'$\mathcal{M}=1$'
    # fields['mach']['lims'] = np.array([1, 2])
    # fields['mach']['cmap'] = 'PiYG'

    for field_name, field in fields.items():
        print(f"Plotting {field_name}")

        fig = matplotlib.figure.Figure()
        axis = fig.add_subplot()

        axis.set_xscale('log')
        axis.set_yscale('log')

        # Axes limits
        axis.set_xlim(1e-1, star.r[-1].to(unyt.r_sun))
        axis.set_ylim(
            companion_masses[0].to(unyt.m_jup),
            companion_masses[-1].to(unyt.m_jup)
        )
        axis.plot(ram_disruption_separations,
                  companion_masses.to(unyt.m_jup), ls='dashed', label='destruction', color='white')
        axis.plot(tidal_disruption_separations,
                  companion_masses.to(unyt.m_jup), ls='dashed', label='destruction')
        axis.plot(destruction_separations,
                  companion_masses.to(unyt.m_jup), label='destruction', ls='dotted', color='black')

        ax2 = axis.twiny()
        ax2.set_xlim((1e-1 * unyt.r_sun).to(unyt.au),
                     star.r[-1].to(unyt.au))
        ax2.set_xscale('log')
        ax2.set_xlabel(
            r'Orbital separation $\ls\unit{\astronomicalunit}\rs$',
            labelpad=10
        )

        # Axes labels
        axis.set_xlabel(r'Orbital separation $\ls R_\odot\rs$')
        axis.set_ylabel(r'Companion mass $\ls M_\text{Jup}\rs$')

        # Axes ticks
        axis.minorticks_on()

        # Plot
        x = star.r.to(unyt.r_sun)[min_idx:]
        y = companion_masses.to(unyt.m_jup)
        z = field['array'].T
        im = axis.pcolormesh(
            x, y, z,
            shading='nearest',
            edgecolors='face',
            cmap=field['cmap'],
            vmin=field['lims'][0], vmax=field['lims'][1],
            rasterized=True
        )

        # eps_rho = 1 contour
        xC, yC = np.meshgrid(x, y)
        contour = axis.contour(
            xC[:, :-500],
            yC[:, :-500],
            epsilon_rho[:-500, :].T,
            [1],
            colors='black', linestyles='dashed'
        )
        axis.clabel(
            contour,
            contour.levels,
            fmt=r"$\varepsilon_\rho=1$",
            manual=[(.3, .3)],
            colors='black',
            fontsize='xx-large',
            use_clabeltext=True,
            inline_spacing=8
        )

        # rp_over_ra = 1 contour
        cont = axis.contour(xC, yC, rp_over_ra.T, [
            1], colors='black', linestyles='dashed')
        axis.clabel(cont, cont.levels, fmt=r"$R_\text{SB} = R_a$", manual=[
            (3, 20)], colors='black', fontsize='xx-large', use_clabeltext=True, inline_spacing=8)

        # a = Rt contour
        # cont2 = ax.contour(xC, yC, a_over_rt.T, [
        #     1], colors='black', linestyles='dotted')
        # ax.clabel(cont2, cont2.levels, fmt=r"$a = R_t$", manual=[
        #     (.5, 2)], colors='black', fontsize='xx-large', use_clabeltext=True, inline_spacing=8)

        # rho v^2 = rho_p vesc^2
        # cont3 = ax.contour(xC, yC, jiaspruit_f.T, [
        #     1], colors='black', linestyles='dotted')
        # ax.clabel(cont3, cont3.levels, fmt=r"$f=1$", manual=[
        #     (.15, 2)], colors='black', fontsize='xx-large', use_clabeltext=True, inline_spacing=8)

        # mach = 1
        # cont4 = ax.contour(xC, yC, mach.T, [1],
        #                    colors='black', linestyles='dashed')
        # ax.clabel(cont4, cont4.levels, fmt=r"$\mathcal{M}=1$", colors='black',
        #           fontsize='xx-large', use_clabeltext=True, inline_spacing=8)

        if field == 'rp_over_ra':
            #    ax.text(20, 20, 'Gravitational regime', fontsize = 'x-large', horizontalalignment = 'center')
            axis.text(3, 30, 'Gravitational regime',
                      fontsize='x-large', horizontalalignment='center')
            axis.text(3, 2, 'Geometrical regime', fontsize='x-large',
                      horizontalalignment='center')

        # Colorbar
        cbar = fig.colorbar(im, pad=0.01)
        cbar.set_label(field['tex_name'], rotation=90)

        fig.savefig(f'plots/parameter_space_{field_name}_{args.radius}.pdf')


main()
