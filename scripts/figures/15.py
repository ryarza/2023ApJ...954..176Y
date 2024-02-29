import numpy as np
import rytools as rt
import matplotlib.pyplot as plt

rt.plot.plottex()
cgs = rt.units.get_cgs()

masses = np.logspace(0, 2, 120) * cgs['MJUP']
r_of_m = rt.planets.mass_radius_relation()
radii = r_of_m(masses / cgs['MJUP']) * cgs['RJUP']
#masses = masses[:2]
data = np.empty_like(masses, dtype = dict)

recompute = False

if recompute:

    for midx, m in enumerate(masses):
        sim_path = '%i/' % midx
        print(midx, m)
        data[midx] = rt.planets.load_orbint_output(sim_path + 'output-total.txt', sim_path + 'scalars.txt', "profiles/rgb_tip/profile1.data", sample = False, nsamples = 100, compute_L = False)

    tau_insp = np.array([data[i]['time'][-1] for i in range(len(masses))])
    eorb_0 = np.array([data[i]['eorb'][0] for i in range(len(masses))])
    dedt_0 = np.array([data[i]['dedt'][0] for i in range(len(masses))])
    vtheta_0 = np.array([data[i]['vtheta'][0] for i in range(len(masses))])
    r_0 = np.array([data[i]['r'][0] for i in range(len(masses))])
    rho_0 = np.array([data[i]['rho'][0] for i in range(len(masses))])
    menc_0 = np.array([data[i]['prof']['itp']['menc'](data[i]['r'])[0] for i in range(len(masses))])
    h_rho_0 = np.array([data[i]['prof']['itp']['h_rho'](data[i]['r'])[0] for i in range(len(masses))])
    rsb = np.array([data[i]['rsb'] for i in range(len(masses))])
    msb = np.array([data[i]['msb'] for i in range(len(masses))])
    ra_0 = 2 * cgs['G'] * msb / pow(vtheta_0, 2)

    print("Radius: ", rsb / cgs['RJUP'])
    edot_0 = - rho_0 * pow(vtheta_0, 3) * np.pi * ( pow(rsb, 2) + pow(ra_0, 2)   )
    v_r_0 = - ( edot_0 / eorb_0 ) * r_0 * ( 1 - 4 * np.pi * pow(r_0, 3) * rho_0 / menc_0 )
    t_low = h_rho_0 / v_r_0
    t_high = r_0 / v_r_0
    print(t_low / cgs['YEAR'])
    print(t_high / cgs['YEAR'])
    print(tau_insp / cgs['YEAR'])
    np.savetxt('tau_insp', tau_insp)

tau_insp = np.loadtxt('tau_insp')

fig, ax = plt.subplots()
ax.plot(masses / cgs['MJUP'], tau_insp / cgs['YEAR'], color = 'black')
ax.plot(masses[-10:] / cgs['MJUP'], tau_insp[-1] * masses[-1] / masses[-10:] / cgs['YEAR'], linestyle = 'dashed', label = r'$\propto M_\text{c}^{-1}$')
ax.plot(masses[:10] / cgs['MJUP'], tau_insp[0] * ( masses[:10] / masses[0] ) * pow(radii[0] / radii[:10], 2) / cgs['YEAR'], linestyle = 'dashed', label = r'$\propto M_\text{c}/R_\text{c}^2$')
ax.set_xlim(1e0, 1e2)
ax.set_ylim(0, 0.12)
ax.tick_params(which = 'both', direction = 'in')
ax.legend(loc = 0)
#ax.set_yscale('log')
#ax.set_xlabel(r'$M_\text{SB}/M_\text{Jup}$')
ax.set_xlabel(r'Companion mass $\ls M_\text{Jup}\rs$')
ax.set_ylabel(r'Orbital decay timescale $\ls\unit{\year}\rs$')
ax.set_xscale('log')
plt.tight_layout()
plt.savefig('test_regime.pdf')
plt.close()
