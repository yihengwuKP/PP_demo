#!/usr/bin/env python

import polymer
import argparse

from openmm import unit

parser = argparse.ArgumentParser()
parser.add_argument('--step-size',   '-ss',  default=0.5,      type=float,
                                            help='step size, default is 0.5 fs')
parser.add_argument('--n-record',   '-nr',  default=1000,   type=int,
                                        help='the total number of frames we want to record')
parser.add_argument('--length',   '-l',  default=10,   type=int,
                                        help='the total number of  carbons in polyethylene.')
parser.add_argument('--force',   '-f',  default=0.0,   type=float,
                                        help='the stretching force on the z axis. Default is 0 kJ/mol/n')
parser.add_argument('--nsteps',     '-ns',  required=True,  type=lambda x: int(float(x)),
                                        help='how many steps should we run the MD for?')
parser.add_argument('--restart',     '-re',  action='store_true', default=False,
                                        help='restart the simulation from checkpnt')

args = parser.parse_args()
print("Initializing...")
fjc = polymer.FJC(args.length)
print("Adding Bonds...")
fjc.add_bond()
if args.force > 0.0:
    fjc.align_z_axis()
    fjc.add_stretching_force(args.force*unit.kilojoules_per_mole/unit.nanometer)
fjc.simulate(n_steps=args.nsteps, n_record=args.n_record, step_size=args.step_size*unit.femtosecond, restart=args.restart)
