"""Compare analytical and numerical orbital decay trajectories."""
import itertools
import os
import shutil
import pathlib
import bookkeeper.nautilus
import rytools.planets
import unyt

os.chdir(pathlib.Path(__file__).parent.parent)

PATH = {}
PATH['GRID'] = pathlib.Path('data/large_sims/inspirals/orbital_trajectories')
PATH['PAR_TEMPLATE'] = PATH['GRID'].parent / "input.txt.template"
PATH['NAUTILUS_ROOT'] = pathlib.Path('/home/rcastroy/src/nautilus')
# PATH['DRAG'] = PATH['NAUTILUS_ROOT'] / 'extras' / 'drag.h5'
PATH['DRAG'] = pathlib.Path('data/drag_uniform_grid.h5')
PATH['MESA_PROFILE_ROOT'] = PATH['NAUTILUS_ROOT'] / 'extras' / 'profiles'


def main():
    """Do main calculation."""
    # for companion_mass in [1, 10, 50, 80] * unyt.mjup:
    for use_drag_coefficients, companion_mass in\
            itertools.product([0, 1], [1, 80] * unyt.mjup):
        companion = rytools.planets.SubstellarBody(companion_mass)
        sim_path = PATH['GRID'] /\
            f"{use_drag_coefficients}" /\
            f"{int(companion_mass.to(unyt.mjup).value)}"
        if sim_path.is_dir():
            shutil.rmtree(sim_path)
        sim_path.mkdir(parents=True, exist_ok=False)
        print(f"Running simulation at {sim_path.resolve()}")
        par = bookkeeper.nautilus.ParameterFile(PATH['PAR_TEMPLATE'])
        par['use_drag_coefficients'] = use_drag_coefficients
        par['drag_data_path'] = str(PATH['DRAG'].resolve())
        par['msb'] = float(companion.mass.in_cgs().value)
        par['rsb'] = float(companion.radius.in_cgs().value)
        if companion_mass == 1 * unyt.m_jup:
            par['mesa_profile_path'] =\
                str((PATH['MESA_PROFILE_ROOT'] / "10rsun.h5").resolve())
        elif companion_mass == 80 * unyt.m_jup:
            par['mesa_profile_path'] =\
                str((PATH['MESA_PROFILE_ROOT'] / "100rsun.h5").resolve())
        par.write(sim_path / "input.txt")
        sim = bookkeeper.nautilus.Simulation(sim_path)
        sim.run(
            [
                str(PATH['NAUTILUS_ROOT'] / 'build' / 'nautilus'),
                sim.par.path.resolve()
            ]
        )


main()
