"""Compare analytical and numerical luminosities."""
import os
import shutil
import pathlib
import bookkeeper.nautilus
import rytools.planets
import unyt

os.chdir(pathlib.Path(__file__).parent.parent)

PATH = {}
PATH['GRID'] = pathlib.Path(
    'data/large_sims/inspirals/compare_energy_deposition_rates')
PATH['PAR_TEMPLATE'] = PATH['GRID'].parent / "input.txt.template"
PATH['NAUTILUS_ROOT'] = pathlib.Path('/home/rcastroy/src/nautilus')
PATH['DRAG'] = pathlib.Path('data') / 'drag_uniform_grid.h5'
PATH['MESA_PROFILES'] =\
    [
        PATH['NAUTILUS_ROOT'] / 'extras' / 'profiles' / '10rsun.h5',
        # PATH['NAUTILUS_ROOT'] / 'extras' / 'profiles' / 'rgb-tip.h5'
]


def main():
    """Do main calculation."""

    assert PATH['DRAG'].is_file(), "Drag file doesn't exist"
    for path in PATH['MESA_PROFILES']:
        assert path.is_file(), "MESA profile doesn't exist"

    for idx_profile, profile in enumerate(PATH['MESA_PROFILES']):
        for companion_mass in [1, 4, 20, 80] * unyt.mjup:
            # for companion_mass in [80] * unyt.mjup:
            companion = rytools.planets.SubstellarBody(companion_mass)
            for use_drag_coefficients in [0, 1]:
                sim_path = PATH['GRID'] /\
                    f"{idx_profile}" /\
                    f"{use_drag_coefficients}" /\
                    f"{int(companion_mass.to(unyt.mjup).value)}"
                if sim_path.is_dir():
                    shutil.rmtree(sim_path)
                sim_path.mkdir(parents=True, exist_ok=False)
                print(f"Running simulation at {sim_path.resolve()}")
                par = bookkeeper.nautilus.ParameterFile(PATH['PAR_TEMPLATE'])
                par['use_drag_coefficients'] = use_drag_coefficients
                par['mesa_profile_path'] = str(profile.resolve())
                par['drag_data_path'] = str(PATH['DRAG'].resolve())
                par['msb'] = float(companion.mass.in_cgs().value)
                par['rsb'] = float(companion.radius.in_cgs().value)
                if companion.mass > 70 * unyt.m_jup:
                    par['output_freq_factor'] = 1.e-3
                else:
                    par['output_freq_factor'] = 1.e-2
                par.write(sim_path / "input.txt")
                sim = bookkeeper.nautilus.Simulation(sim_path)
                sim.run(
                    [
                        str(PATH['NAUTILUS_ROOT'] / 'build' / 'nautilus'),
                        str(sim.par.path.resolve())
                    ]
                )


main()
