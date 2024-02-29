"""Plots convergence of ram pressure drag with resolution at boundary."""
import pathlib
import h5py
import matplotlib
import matplotlib.figure
import matplotlib.style
import numpy as np
import rytools as rt
from rytools import cewt

rt.set_work_dir(pathlib.Path(__file__))

matplotlib.style.use('~/.config/matplotlib/ricardo.mplstyle')

paths = [pathlib.Path(f'large_sims/test-resolution-at-boundary/max_{i}')
         for i in [2, 3, 4, 5]]
fac = np.array([pow(2, i) for i in [1, 2, 3, 4]])
dx_min = (10 / 128) / fac
n_cells = 0.6 / dx_min
p_drag = np.empty_like(paths, dtype=float)
g_drag = np.empty_like(paths, dtype=float)

for idx, path in enumerate(paths):
    p_drag[idx] = cewt.load_cp(path / 'ram-pressure.h5')
    g_drag[idx] = cewt.load_cg_av(path, 1.6, time_av_low=20)
    print(f"Pressure drag for {path}: {p_drag[idx]:.5e}")
    print(f"Grav drag for {path}: {g_drag[idx]:.5e}")

fig = matplotlib.figure.Figure()
ax = fig.add_subplot()
eps = np.abs(p_drag / p_drag[-1] - 1)
ax.set_xlabel(r'$R_\text{SB}/\Delta x_\text{min}$')
ax.set_ylabel(r'Relative error')
ax.set_yscale('log')
ax.set_ylim(1e-3, 1e0)
print(eps)
ax.plot(n_cells[:-1], eps[:-1], label='Ram pressure')
# eps = np.abs(g_drag / g_drag[-1] - 1)
# ax.plot(n_cells[:-1], eps[:-1], label='Gravitational')
# print(eps)
print(n_cells)
fig.savefig('plots/resolution_convergence.png')

prop_cycle = matplotlib.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']

fig = matplotlib.figure.Figure()
ax = fig.add_subplot()

for idx, path in enumerate(paths):
    data = h5py.File(path / "ram-pressure.h5", 'r')
    rp_over_ra = data.attrs['r']
    assert rp_over_ra > 0
    assert (np.diff(data['time'][()]) > 0).all(), "Time must increase"

    force = data['ram-pressure-drag'][:, 0]

    # Density and speed are 1
    mask = data['time'][()] <= 20.01
    drag_coefficient = force / np.pi / pow(rp_over_ra, 2)
    ax.plot(data['time'][()][mask] / 0.6, drag_coefficient[mask],
            label=f"{round(n_cells[idx])}"r" cells per $R_\text{SB}$",
            color=colors[idx])

    grav_data = cewt.load_cg(path, 1.6)
#    ax.plot(grav_data["time"], grav_data["drag_coefficient"], ls='dashed')


def xtox2(x):
    """Convert to secondary axis."""
    return x * 0.6


def x2tox(x):
    """Convert to primary axis."""
    return x / 0.6


secax = ax.secondary_xaxis('top', functions=(xtox2, x2tox))
secax.set_xlabel(r'Time $\ls R_a/v_\infty\rs$')
# ax.set_ylim(-0.5, 1.5)
ax.set_ylim(0, 1.5)
ax.set_xlim(0, 20 / 0.6)
ax.autoscale(enable=False)
ax.fill_between([15 / 0.6, data["time"][()][mask][-1] / 0.6],
                [ax.get_ylim()[0]] * 2, [ax.get_ylim()[1]] * 2,
                color='gray', alpha=0.3)
ax.text(17.5 / 0.6, 0.4, "Time-averaging\nwindow",
        ha='center', va='center', fontsize='large')
ax.set_xlabel(r"Time $\ls R_\text{SB}/v_\infty\rs$")
ax.set_ylabel(r"Ram pressure drag coefficient")
# ax.legend(ncol=2, loc='lower right')
ax.legend(loc=0)
fig.savefig('plots/convergence_drag_with_time.pdf')
