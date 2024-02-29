import argparse
import os
import pathlib
import shutil
import warnings

import h5py
import mesa_reader as mr
import numpy as np
import numpy.typing as npt
import scipy
import rytools.planets
import rytools.units
import rytools as rt
from mpi4py import MPI
import bookkeeper.nautilus
import unyt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
MPI_SIZE = comm.Get_size()

os.chdir(pathlib.Path(__file__).parent.parent)


def function_to_optimize(companion_mass: float, profile: rt.stars.MesaProfile):
    # Create directory for simulation
    companion = rt.planets.SubstellarBody(companion_mass * unyt.m_jup)
    sim_path = pathlib.Path(
        f"data/outcomes_grid/{profile.r[-1].to(unyt.r_sun).value:.5e}"
        f"/{companion.mass.to(unyt.m_jup).value:.5e}"
    )
    if sim_path.is_dir():
        shutil.rmtree(sim_path)
    sim_path.mkdir(parents=True)
    # Copy profile to simulation directory and convert it to readable format
    shutil.copyfile(profile.path, sim_path / "profile.data")
    sim_profile = rt.stars.MesaProfile(sim_path / "profile.data")
    sim_profile.save_to_hdf5(sim_path / "profile.h5", include_hash=True)

    # Copy and modify input.txt template
    shutil.copyfile(
        'data/outcomes_grid/input.txt.template',
        sim_path / "input.txt"
    )

    sim = bookkeeper.nautilus.Simulation(sim_path)
    sim.par['msb'] = float(companion.mass.to(unyt.g).value)
    sim.par['rsb'] = float(companion.radius.to(unyt.cm).value)
    sim.par.write()

    # print("Running sim...")
    sim.run(
        run_command=['/home/rcastroy/src/nautilus/build/nautilus', 'input.txt']
    )

    sim_output = rt.planets.EngulfmentIntegratorSimulation(sim_path)
    # print(sim_output.companion.mass.to(unyt.m_jup))
    # print(sim_output.prof.r[-1].to(unyt.r_sun))
    if sim_output.envelope_ejection:
        return 0.5
    return -0.5


def function_to_optimize_conv_region(companion_mass: float, profile: rt.stars.MesaProfile):
    # Create directory for simulation
    companion = rt.planets.SubstellarBody(companion_mass * unyt.m_jup)
    sim_path = pathlib.Path(
        f"data/outcomes_grid/{profile.r[-1].to(unyt.r_sun).value:.6e}"
        f"/{companion.mass.to(unyt.m_jup).value:.6e}"
    )
    if sim_path.is_dir():
        shutil.rmtree(sim_path)
    sim_path.mkdir(parents=True)
    # Copy profile to simulation directory and convert it to readable format
    shutil.copyfile(profile.path, sim_path / "profile.data")
    sim_profile = rt.stars.MesaProfile(sim_path / "profile.data")
    sim_profile.save_to_hdf5(sim_path / "profile.h5", include_hash=True)

    # Copy and modify input.txt template
    shutil.copyfile(
        'data/outcomes_grid/input.txt.template',
        sim_path / "input.txt"
    )

    sim = bookkeeper.nautilus.Simulation(sim_path)
    sim.par['msb'] = float(companion.mass.to(unyt.g).value)
    sim.par['rsb'] = float(companion.radius.to(unyt.cm).value)
    sim.par.write()

    # print("Running sim...")
    sim.run(
        run_command=['/home/rcastroy/src/nautilus/build/nautilus', 'input.txt']
    )

    sim_output = rt.planets.EngulfmentIntegratorSimulation(sim_path)
    destroyed_separation = sim_output.r[sim_output.destroyed_idx]
    base_of_convective_zone = sim_output.prof.r[sim_output.prof.base_of_convective_zone_index]
    # print(f"star radius {profile.r[-1]}")
    # print(f"{base_of_convective_zone}")
    # exit(0)
    if destroyed_separation < base_of_convective_zone:
        # print(f"Mass {companion.mass} destroyed below base")
        return 0.5
    # print(f"Mass {companion.mass} destroyed above base")
    return -0.5


parser = argparse.ArgumentParser()
parser.add_argument(
    "-f", "--folder", help="Folder with MESA output files",
    required=False,
    default=pathlib.Path
    (
        """/data/groups/ramirez-ruiz/rcastroy/mesa/grid_sun_like/"""
        """1.00000e+00/LOGS_to_start_he_core_flash"""
    )
)
args = parser.parse_args()

# Planet masses (cgs)
companion_masses = np.geomspace(1, 100, 200) * unyt.m_jup
companions = [rytools.planets.SubstellarBody(mass)
              for mass in companion_masses]


def get_profile_numbers(n_profiles: int) -> npt.NDArray:
    logdir = mr.MesaLogDir(str(args.folder.resolve()))

    assert len(logdir.profile_numbers) == len(logdir.history.R)

    sim_radii = logdir.history.R

    sim_radii_increasing = rt.filter_array_increasing(sim_radii)

    desired_radii = np.geomspace(
        max(2, logdir.history.R[0]),
        np.amax(logdir.history.R),
        n_profiles - 1
    )
    desired_radii = np.append(desired_radii, 11.43)
    desired_radii = np.sort(desired_radii)

    if rank == 0:
        print(f"  Analyzing {n_profiles} profiles between"
              f"R = {desired_radii[0]:.3f}"
              f"RSun and {desired_radii[-1]:.3f} RSun")

    profile_numbers = np.empty_like(desired_radii, dtype=int)
    for idx, desired_radius in enumerate(desired_radii):
        prof_idx, val1 = rt.nearest(sim_radii_increasing, desired_radius)
        if not np.isclose(desired_radius, val1, rtol=1.e-2) and rank == 0:
            warnings.warn(
                f"Desired radius {desired_radius:.10e}, but closest is"
                f"profile {logdir.profile_numbers[prof_idx]}"
                f"with radius {val1:.10e}")
        prof_idx, val2 = rt.nearest(sim_radii, val1)
        assert np.isclose(val1, val2, rtol=1.e-12)
        profile_numbers[idx] = logdir.profile_numbers[prof_idx]

    return profile_numbers


def main() -> None:
    n_profiles = 30
    profile_numbers = get_profile_numbers(n_profiles)
    h = mr.MesaData(str(args.folder.resolve() / 'history.data'))
    # profile_idxs = rt.filter_array_increasing_indices(l.history.R)
    # n_profiles = len(profile_idxs)
    stellar_radii = h.R[profile_numbers] * unyt.r_sun
    stellar_ages = h.star_age[profile_numbers] * unyt.yr

    assert (np.diff(stellar_radii) > 0).all(), np.amin(np.diff(stellar_radii))
    assert (np.diff(stellar_ages) > 0).all()

    if rank == 0:
        print(f"Total number of profiles: {n_profiles}")
        minimum_ejection_masses = np.empty(n_profiles, dtype=float)
        minimum_survive_conv_masses = np.empty(n_profiles, dtype=float)
        # final_separations = np.empty((n_profiles, len(companions)))
        # conv_zone_start_separations = np.empty(n_profiles)

    NSWEEPS = int(np.ceil(n_profiles / MPI_SIZE))
    for sweep in range(NSWEEPS):
        i = MPI_SIZE * sweep + rank
        if i < n_profiles:
            star = rt.stars.MesaProfile(
                args.folder / f'profile{(profile_numbers[i])}.data')

            # Envelope ejection line
            if function_to_optimize(100, star) == 0.5:
                minimum_ejection_mass = scipy.optimize.brentq(
                    function_to_optimize,
                    10, 100,
                    args=(star),
                    rtol=1.e-2,
                    maxiter=100
                )
            else:
                minimum_ejection_mass = np.nan

            # Destroyed below convective zone line
            if function_to_optimize_conv_region(0.1, star) == -0.5:
                minimum_survive_conv_mass = scipy.optimize.brentq(
                    function_to_optimize_conv_region,
                    0.1, 20,
                    args=(star),
                    rtol=1.e-2,
                    maxiter=100
                )
            else:
                minimum_survive_conv_mass = np.nan

            print(
                f"Radius: {star.r[-1]}. Min eject mass: {minimum_ejection_mass}")
            print(
                f"Radius: {star.r[-1]}. Min mass to survive"
                f" below conv region: {minimum_survive_conv_mass}")

            # Gather results in rank 0
            if rank != 0:
                comm.send(minimum_ejection_mass, dest=0, tag=0)
                comm.send(minimum_survive_conv_mass, dest=0, tag=1)
                # comm.send(final_separation, dest=0, tag=1)
                # comm.send(prof['r'][prof['conv_zone_edge']], dest=0, tag=2)
            else:
                if sweep < NSWEEPS - 1:
                    number_of_senders = MPI_SIZE
                else:
                    number_of_senders = MPI_SIZE - \
                        (NSWEEPS * MPI_SIZE - n_profiles)
                for sender in range(1, number_of_senders):
                    minimum_ejection_masses[MPI_SIZE * sweep + sender] =\
                        comm.recv(source=sender, tag=0)
                    minimum_survive_conv_masses[MPI_SIZE * sweep + sender] =\
                        comm.recv(source=sender, tag=1)
                print("Rank 0 finished receiving mass indices for sweep %i out of %i" % (
                    sweep + 1, NSWEEPS))
                minimum_ejection_masses[MPI_SIZE * sweep] =\
                    minimum_ejection_mass
                minimum_survive_conv_masses[MPI_SIZE * sweep] =\
                    minimum_survive_conv_mass
                # final_separations[MPI_SIZE * sweep] = final_separation
                # conv_zone_start_separations[MPI_SIZE *
                # sweep] = prof['r'][prof['conv_zone_edge']]

    if rank == 0:
        f = h5py.File('data/outcomes.h5', 'w')
        f.create_dataset(
            'minimum_survive_conv_masses',
            data=minimum_survive_conv_masses
        )
        f.create_dataset(
            'minimum_ejection_masses',
            data=minimum_ejection_masses
        )
        f.create_dataset('stellar_radii', data=stellar_radii.to(unyt.r_sun))
        f.close()


main()
