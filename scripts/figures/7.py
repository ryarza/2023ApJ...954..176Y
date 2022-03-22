import numpy as np
import matplotlib.pyplot as plt
import rytools as rt
import argparse
import h5py

rt.plot.plottex()

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputdata", help = "Path to drag coefficients", required = True)
args = parser.parse_args()

# Values of epsilon_rho and Rp/Ra used in the grid
q = np.logspace(-3, -1, 4)[2]
eps_rhos = np.logspace(-1, 0, 4)
rp_over_ras = np.logspace(np.log10(0.3), 0, 7)

def formats(number):
    if number == int(number):
        return "%i" % number
    else:
        return "%.1f" % number

d = rt.h5_to_dict(args.inputdata)

markers = ['o', 'x', 's', '^']

# Plot gravitational drag
fig, ax = plt.subplots(constrained_layout = True)
for idx_eps_rho, eps_rho in enumerate(eps_rhos):
    ax.plot(rp_over_ras, d['cg'][2][idx_eps_rho], label = r'$\varepsilon_\rho = %.2f$' % eps_rho, marker = markers[idx_eps_rho])
ax.tick_params(which = 'both', direction = 'in')
ax.set_xlabel(r'$R_\text{SB} / R_a$')
ax.set_ylabel(r'$C_g$')
ax.legend(loc = 0)
plt.savefig('7_right.pdf', bbox_inches = 'tight')
plt.close()

# Plot geometrical drag
fig, ax = plt.subplots(constrained_layout = True)
for idx_eps_rho, eps_rho in enumerate(eps_rhos):
    ax.plot(rp_over_ras, d['cp'][2][idx_eps_rho] * pow(rp_over_ras, 2), label = r'$\varepsilon_\rho = %.2f$' % eps_rho, marker = markers[idx_eps_rho])
ax.tick_params(which = 'both', direction = 'in')
ax.scatter(rp_over_ras[-1], np.array([0.25]), label = 'Smooth sphere', color = 'black', marker = 'd')
ax.set_xlabel(r'$R_\text{SB} / R_a$')
ax.set_ylabel(r'$C_p\lp R_\text{SB} / R_a\rp^2$')
ax.legend(loc = 0)
plt.savefig('7_left.pdf', bbox_inches = 'tight')
plt.close()
