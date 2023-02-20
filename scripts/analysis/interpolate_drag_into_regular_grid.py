"""Interpolate drag coefficients onto a regular grid in parameter space."""
import itertools
import os
import pathlib
import h5py
import numpy as np
import pandas as pd
# import numpy.typing as npt
# import scipy
import rytools as rt
import scipy

os.chdir(pathlib.Path(__file__).parent)


def regular_grid_points(dataset: dict) -> tuple:
    """Regular points onto which to interpolate the drag coefficients."""
    df = pd.DataFrame()
    npoints = 100
    mass_ratio = np.geomspace(
        dataset["mass_ratio"].min(),
        dataset["mass_ratio"].max(),
        npoints
    )

    eps_rho = np.geomspace(
        dataset["eps_rho"].min(),
        dataset["eps_rho"].max(),
        npoints
    )

    rp_over_ra = np.geomspace(
        dataset["rp_over_ra"].min(),
        dataset["rp_over_ra"].max(),
        npoints
    )

    product = np.log10(np.array(
        list(itertools.product(mass_ratio, eps_rho, rp_over_ra))))

    df['log_mass_ratio'] = product[:, 0]
    df['log_eps_rho'] = product[:, 1]
    df['log_rp_over_ra'] = product[:, 2]

    unique = {}
    unique['log_mass_ratio'] = np.log10(mass_ratio)
    unique['log_eps_rho'] = np.log10(eps_rho)
    unique['log_rp_over_ra'] = np.log10(rp_over_ra)

    return df, product, unique


def main() -> None:
    """Do main calculation."""
    save_path = pathlib.Path('data/drag_uniform_grid.h5')
    data_path = pathlib.Path('data/drag.h5')
    print(f"Loading drag coefficients from {data_path.resolve()}")
    data = rt.h5_to_dict(pathlib.Path(data_path))

    c_ram_eff = data['c_ram'] * np.minimum(1, pow(data["rp_over_ra"], 2))
    c_grav_eff = data['c_grav'] * np.minimum(1, pow(data["rp_over_ra"], -2))

    assert (c_ram_eff + c_grav_eff > 0).all(), "Experiencing negative drag!"
    # print(np.amin(c_ram_eff))
    # print(np.amin(c_grav_eff))
    # for idx in range(len(c_ram_eff)):
    #     if (c_ram_eff[idx] + c_grav_eff[idx]) < 0:
    #         # if c_ram_eff[idx] < -0.01:
    #         print(data["mass_ratio"][idx], data["eps_rho"]
    #               [idx], data["rp_over_ra"][idx], data['c_ram'][idx], c_ram_eff[idx])
    #         exit(0)

    unstructured_points = np.log10(
        np.array(
            [data["mass_ratio"], data["eps_rho"], data["rp_over_ra"]]
        ).T
    )

    regular_points, product, unique = regular_grid_points(data)

    # Linear interpolation for points inside domain
    print("Interpolating regular points...")
    regular_points['c_ram'] = scipy.interpolate.griddata(
        unstructured_points,
        data["c_ram"],
        product,
        rescale=True
    )

    regular_points['c_grav'] = scipy.interpolate.griddata(
        unstructured_points,
        data["c_grav"],
        product,
        rescale=True
    )

    # Nearest interpolation for points outside simulation domain
    print("Interpolating NaN points...")
    mask = np.isnan(regular_points['c_grav'])
    assert (mask == np.isnan(regular_points['c_ram'])).all()

    nan_points_itp_format = product[mask]

    regular_points.loc[mask, 'c_grav'] =\
        scipy.interpolate.griddata(
        unstructured_points,
        data["c_grav"],
        nan_points_itp_format,
        method='nearest',
        rescale=True
    )

    regular_points.loc[mask, 'c_ram'] =\
        scipy.interpolate.griddata(
        unstructured_points,
        data["c_ram"],
        nan_points_itp_format,
        method='nearest',
        rescale=True
    )

    # As a result of interpolation, some points will have negative drag.
    # Use nearest interpolation for those points too

    # c_ram_eff = regular_points['c_ram'] * \
    #     np.minimum(1, pow(pow(10, regular_points["log_rp_over_ra"]), 2))
    # c_grav_eff = regular_points['c_grav'] * \
    #     np.minimum(1, pow(pow(10, regular_points["log_rp_over_ra"]), -2))

    # mask = np.where(c_ram_eff + c_grav_eff < 0)[0]

    # regular_points.loc[mask, 'c_grav'] =\
    #     scipy.interpolate.griddata(
    #     unstructured_points,
    #     data["c_grav"],
    #     product[mask],
    #     method='nearest',
    #     rescale=True
    # )

    # regular_points.loc[mask, 'c_ram'] =\
    #     scipy.interpolate.griddata(
    #     unstructured_points,
    #     data["c_ram"],
    #     product[mask],
    #     method='nearest',
    #     rescale=True
    # )

    # for i in range(len(regular_points)):
    #     if np.isclose(regular_points['log_mass_ratio'][i], np.log10(4.e-3), rtol=1e-2):
    #         if np.isclose(regular_points['log_eps_rho'][i], -1, rtol=1e-2):
    #             print(regular_points['log_mass_ratio'][i], regular_points['log_eps_rho'][i],
    #                   regular_points['log_rp_over_ra'][i], regular_points['c_ram'][i], regular_points['c_grav'][i])

    assert not np.isnan(regular_points['c_ram']).any()
    assert not np.isnan(regular_points['c_grav']).any()

    # c_ram_eff = regular_points['c_ram'] * \
    #     np.minimum(1, pow(pow(10, regular_points["log_rp_over_ra"]), 2))
    # c_grav_eff = regular_points['c_grav'] * \
    #     np.minimum(1, pow(pow(10, regular_points["log_rp_over_ra"]), -2))

    # for i in range(len(c_ram_eff)):
    #     if not (c_ram_eff[i] + c_grav_eff[i] > 0):
    #         print(
    #             pow(10, regular_points["log_mass_ratio"][i]),
    #             pow(10, regular_points["log_eps_rho"][i]),
    #             pow(10, regular_points["log_rp_over_ra"][i]),
    #             f"c_ram = {regular_points['c_ram'][i]} ",
    #             f"c_grav = {regular_points['c_grav'][i]} ",
    #             c_ram_eff[i] + c_grav_eff[i]
    #         )
    #         assert regular_points['c_ram'][i] in data['c_ram']
    # print("Experiencing negative drag!")
    # print(np.amin(c_ram_eff + c_grav_eff))

    print("Saving drag coefficients on regular grid to "
          f"{save_path.resolve()}")
    with h5py.File(save_path, 'w') as dset:

        dset.attrs['time_av_low'] = data['time_av_low']
        dset.attrs['integration_radius'] = data['integration_radius']
        dset.create_dataset(
            'npoints',
            data=np.array(
                [len(unique[i])
                 for i in ['log_mass_ratio', 'log_eps_rho', 'log_rp_over_ra']],
                dtype=int
            )
        )

        for itp_var, itp_vals in unique.items():
            dset.create_dataset(itp_var, data=itp_vals)

        c_grav_reshaped = np.array(regular_points['c_grav']).reshape(
            dset['npoints'][()]
        )

        # for i, logq in enumerate(unique['log_mass_ratio']):
        #     for j, logeps in enumerate(unique['log_eps_rho']):
        #         for k, logr in enumerate(unique['log_rp_over_ra']):
        #             point = (logq, logeps, logr)
        #             val = scipy.interpolate.griddata(
        #                 unstructured_points,
        #                 data["c_grav"],
        #                 point
        #             )
        #             if not np.isnan(val):
        #                 if not np.isclose(val, c_grav_reshaped[i, j, k], rtol=1e-1):
        #                     print(i, j, k, val, c_grav_reshaped[i, j, k])

        c_ram_reshaped = np.array(regular_points['c_ram']).reshape(
            dset['npoints'][()]
        )

        dset.create_dataset('c_grav', data=c_grav_reshaped)
        dset.create_dataset('c_ram', data=c_ram_reshaped)


main()
