#!/usr/bin/python3

import f90nml
import numpy as np
import h5py
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

plt.rcParams["font.family"] = "serif"
plt.rcParams.update({'font.size': 9})


def Readfreq(filename) :

    print("....... Reading FC2 file: {0} ........".format(filename))
    h5f = h5py.File(filename, mode='r')

    dsetname = 'freq_dat'
    freq_dat = h5f[dsetname][:]
    h5f.close()

    return freq_dat

def Plot(freq_dat, bais_atm, num_pnt, prefix, T) :

    print("Plotting ...")
    fig = plt.figure()
    ax = fig.add_subplot()
    clr = 'black'

    eps = 1.0E-3
    x = freq_dat[:, 0]
    xmin = np.min(x[:]) - eps
    xmax = np.max(x[:]) + eps

    dof = 3*basis_atm
    for k in range(1, dof+1) :

        y = freq_dat[:, k]

        ax.plot(x, y, '-', linewidth=1.0, color=clr)

    ymin = np.min(freq_dat[:, 1:]) - eps
    ymax = np.max(freq_dat[:, 1:]) + 0.1

    strt = 0
    x_tick = []
    for nn in range(1, len(num_pnt)) :

        end = strt + num_pnt[nn]
        freq_part = freq_dat[strt:end, :]
        strt += num_pnt[nn]
        x_tick.append(freq_dat[end-1, 0])

        ax.vlines(x=freq_dat[end-1, 0], ymin=ymin, ymax=ymax, colors='gray', ls='--', lw=1.0)

    ax.set_xticks(x_tick)
    ax.set_xticklabels([])

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.yaxis.set_minor_locator(MultipleLocator(20))
    ax.tick_params(axis="y",direction="in", length=5, width=1)
    ax.tick_params(which='minor', direction="in", length=3, width=0.8)
    ax.tick_params(axis="y", which='minor', direction="in", length=3, width=0.8)
    ax.yaxis.set_ticks_position('both')
    #ax.yaxis.tick_right()
    
    ax.set_xlabel('q', fontsize=12)
    ax.set_ylabel(r'$\nu \; (THz)$', fontsize=12)
    for label in ax.get_xticklabels() :
        label.set_fontsize(12)
    
    fname = prefix+'_phonon_disp_'+str(T)+'K.pdf'
    plt.savefig(fname, format="pdf", bbox_inches="tight")
    plt.show()



if __name__ == "__main__" :

    parser = argparse.ArgumentParser()

    parser.add_argument("-inp", type=str, default='input.nml', \
            help="Enter the name of input file")
    parser.add_argument("-T", type=float, default=300, \
            help="Temperature in Kelvin")

    args = parser.parse_args()
    fname = args.inp
    T = args.T

    data_parser = f90nml.Parser()
    nml_data = data_parser.read(fname)

    prefix = nml_data['FCInfo']['prefix']
    basis_atm = nml_data['CrystalInfo']['Nbasis']
    num_pnt = np.array(nml_data['HighSymPath']['q_high_sym'], dtype=np.intc)[:, -1]

    filename = 'dispersion_'+prefix+str(T)+'K.h5'
    freq_dat = Readfreq(filename)
    
    Plot(freq_dat, basis_atm, num_pnt, prefix, T)

