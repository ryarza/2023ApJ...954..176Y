"""Interpolate drag coefficients onto a regular grid in parameter space."""
import os
import pathlib
import h5py
from rytools import cewt

os.chdir(pathlib.Path(__file__).parent)


def main() -> None:
    """Do main calculation."""
    time_av_low = 15
    integration_radius = 1.6
    save_path = pathlib.Path('data/drag.h5')
    grid = cewt.SimulationGrid(pathlib.Path('data/large_sims/parameter_space'))

    # Choose only simulations that have ram pressure computed
    paths = [sim.path for sim in grid if
             (sim.path / "ram_pressure.h5").is_file() and sim.complete]

    grid = cewt.SimulationGrid(paths)

    print(f"Computing drag for {len(grid)} simulations")

    c_grav = grid.c_grav_av(time_av_low, integration_radius)
    c_ram = grid.c_ram_av(time_av_low)

    print(f"Saving drag coefficients to file {save_path}")

    dset = h5py.File(save_path, 'w')

    dset.attrs['time_av_low'] = time_av_low
    dset.attrs['integration_radius'] = integration_radius

    dset.create_dataset('mass_ratio', data=grid['sim_q'])
    dset.create_dataset('eps_rho', data=grid.eps_rho)
    dset.create_dataset('rp_over_ra', data=grid['sim_radius'])
    dset.create_dataset('c_ram', data=c_ram)
    dset.create_dataset('c_grav', data=c_grav)
    dset.close()

    # df = pd.DataFrame()
    # df.attrs['time_av_low'] = time_av_low
    # df.attrs['integration_radius'] = integration_radius

    # df['mass_ratio'] = grid['sim_q']
    # df['eps_rho'] = grid.eps_rho
    # df['rp_over_ra'] = grid['sim_radius']
    # df['c_ram'] = c_ram
    # df['c_grav'] = c_grav

    # # df.to_csv(save_path)
    # df.to_hdf(save_path, 'main', mode='w')


main()
