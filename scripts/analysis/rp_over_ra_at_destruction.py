import argparse
import pathlib
import numpy as np
import mesa_reader as mr
import rytools
import rytools.planets
from mpi4py import MPI
import h5py
import unyt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
MPI_SIZE = comm.Get_size()

parser = argparse.ArgumentParser()
parser.add_argument(
    "-f", "--folder",
    help="Folder with MESA output files",
    type=pathlib.Path,
    default=pathlib.Path(
        '/data/groups/ramirez-ruiz/rcastroy/mesa/'
        'grid_sun_like/1.00000e+00/LOGS_to_start_he_core_flash/'
    )
)
parser.add_argument(
    "-o", "--outputfile",
    help="Where to store output",
    type=pathlib.Path,
    default=pathlib.Path('data/rp_over_ra_at_destruction.h5')
)
# parser.add_argument("-l", "--logdefault",
#                     help="If folder isn't available, try LOGS instead", required=True)
args = parser.parse_args()
# args.logdefault = args.logdefault == "True"

if rank == 0:
    print(f"Working on folder {args.folder.resolve()}")
if not args.folder.is_dir():
    if rank == 0:
        print(f"  Folder {args.folder.resolve()} doesn't exist")
    # if args.logdefault:
    #     args.folder = args.folder.parent / 'LOGS'
    #     if rank == 0:
    #         print(f"  Working on folder {args.folder.resolve()} instead")
    # else:
    exit(0)

mps = np.linspace(10, 80, 200) * unyt.m_jup
sbs = [rytools.planets.SubstellarBody(mass) for mass in mps]

# Desired radii
ld = mr.MesaLogDir(args.folder)

assert len(ld.profile_numbers) == len(ld.history.R)

sim_radii = ld.history.R

sim_radii_increasing = np.array([sim_radii[0]], dtype=float)
for r in sim_radii:
    if r > sim_radii_increasing[-1]:
        sim_radii_increasing = np.append(sim_radii_increasing, r)

# desired_radii = sim_radii_increasing
n_desired = 200
desired_radii = np.geomspace(
    max(10, ld.history.R[0]),
    np.amax(ld.history.R),
    n_desired
)

if rank == 0:
    print(f"  Analyzing {n_desired} profiles between"
          f"R = {desired_radii[0]:.3f} RSun and {desired_radii[-1]:.3f} RSun")

profile_numbers = np.empty_like(desired_radii, dtype=int)
for idx, desired_radius in enumerate(desired_radii):
    prof_idx, val1 = rytools.nearest(sim_radii_increasing, desired_radius)
    if not np.isclose(desired_radius, val1, rtol=1.e-2) and rank == 0:
        print("Warning: desired radius %.10e, but closest is profile %i with radius %.10e" % (
            desired_radius, ld.profile_numbers[prof_idx], val1))
    prof_idx, val2 = rytools.nearest(sim_radii, val1)
    assert np.isclose(val1, val2, rtol=1.e-12)
    profile_numbers[idx] = ld.profile_numbers[prof_idx]

n_profiles = len(desired_radii)

NSWEEPS = int(np.ceil(n_profiles / MPI_SIZE))

if rank == 0:
    print(f"  Number of sweeps: {NSWEEPS}")
    rp_over_ras = np.empty((len(desired_radii), len(mps)))
    minimum_ejection_masses_alpha = np.empty(n_profiles)

rp_over_ras_local = np.empty_like(mps)

for sweep in range(NSWEEPS):
    i = MPI_SIZE * sweep + rank
    if i < n_profiles:
        star = rytools.stars.MesaProfile(
            args.folder / f'profile{profile_numbers[i]}.data'
        )

        for idx, sb in enumerate(sbs):
            ce = rytools.planets.AnalyticalCommonEnvelope(star, sb)
            rp_over_ras_local[idx] = ce.rsb_over_ra[ce.companion_destruction_idx]

    #     Gather results in rank 0
        if rank != 0:
            comm.send(rp_over_ras_local, dest=0, tag=0)
        else:
            if sweep < NSWEEPS - 1:
                number_of_senders = MPI_SIZE
            else:
                number_of_senders = MPI_SIZE - \
                    (NSWEEPS * MPI_SIZE - n_profiles)
            for sender in range(1, number_of_senders):
                rp_over_ras[MPI_SIZE * sweep +
                            sender] = comm.recv(source=sender, tag=0)
            print("  Rank 0 finished receiving mass indices for sweep %i out of %i" % (
                sweep, NSWEEPS - 1))
            rp_over_ras[MPI_SIZE * sweep] = rp_over_ras_local

if rank == 0:
    f = h5py.File(args.outputfile, 'w')
    f.create_dataset('planet_mass', data=mps.to(unyt.m_jup))
    f.create_dataset('stellar_radius', data=desired_radii)
    f.create_dataset('rp_over_ra', data=rp_over_ras)
    f.close()
    print(f"Saved results to {args.outputfile}")
