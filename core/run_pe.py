#!/usr/bin/env python

import pe
import argparse

from openmm import unit

parser = argparse.ArgumentParser()
parser.add_argument('--step-size',   '-ss',  default=0.5,      type=float,
                                            help='step size, default is 0.5 fs')
parser.add_argument('--n-record',   '-nr',  default=1000,   type=int,
                                        help='the total number of frames we want to record')
parser.add_argument('--length',   '-l',  default=10,   type=int,
                                        help='the total number of  carbons in polyethylene.')
parser.add_argument('--nsteps',     '-ns',  required=True,  type=lambda x: int(float(x)),
                                        help='how many steps should we run the MD for?')

args = parser.parse_args()
pe = pe.PE(args.nsteps)
pe.add_bond()
pe.add_angle()
pe.add_dihedral()
pe.add_LJ()
pe.simulate(args.nsteps, args.n_record, step_size=args.step_size*unit.femtosecond)
