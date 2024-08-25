import openmm as mm
from openmm import app, unit
import networkx as nx
import numpy as np
from itertools import combinations, product


mass = {"C": 12.011, "H": 1.008}


class PE():
    name2fg = {'bond': 0, 'angle': 1, 'dihedral': 2, 'LJ': 3}
    fg2name = {0: 'bond', 1: 'angle', 2: 'dihedral', 3: 'LJ'}

    def __init__(self, length=10):
        self.init_config(N_carbons=length)
        self.init_topology(N_carbons=length)
        #self.write_psf()
        self.write_init_config()
        self.init_system()
        self.fgs = []

    def init_config(self, N_carbons=10):

        graph = nx.Graph()

        theta = 109.5 / 180 * np.pi / 2
        b_CH = 1.09
        b_CC = 1.526
        pos = []
        pos.append([0, -b_CH*np.sin(theta), b_CH*np.cos(theta)])   # H on initial C
        pos.append([0, 0, 0])  # first C
        pos.append([b_CH*np.sin(theta), 0, -b_CH*np.cos(theta)])   # H on initial C
        pos.append([-b_CH*np.sin(theta), 0, -b_CH*np.cos(theta)])  # H on initial C
        graph.add_node(0, name='H')
        graph.add_node(1, name='C')
        graph.add_node(2, name='H')
        graph.add_node(3, name='H')
        graph.add_edge(1, 0)
        graph.add_edge(1, 2)
        graph.add_edge(1, 3)

        index = 4
        for i in range(1, N_carbons):
            if i % 2 == 1:
                pos.append([0, i*b_CC*np.sin(theta), b_CC*np.cos(theta)])
                pos.append([b_CH*np.sin(theta), i*b_CC*np.sin(theta), (b_CC+b_CH)*np.cos(theta)])
                pos.append([-b_CH*np.sin(theta), i*b_CC*np.sin(theta), (b_CC+b_CH)*np.cos(theta)])
                graph.add_node(index, name='C')
                graph.add_node(index+1, name='H')
                graph.add_node(index+2, name='H')
                graph.add_edge(index-3, index)
                graph.add_edge(index, index+1)
                graph.add_edge(index, index+2)
                index += 3
                if i == N_carbons-1:
                    pos.append([0, i*b_CC*np.sin(theta)+b_CH*np.sin(theta), (b_CC-b_CH)*np.cos(theta)])
                    graph.add_node(index, name='H')
                    graph.add_edge(index-3, index)

            else:
                pos.append([0, i*b_CC*np.sin(theta), 0])
                pos.append([b_CH*np.sin(theta), i*b_CC*np.sin(theta), -b_CH*np.cos(theta)])
                pos.append([-b_CH*np.sin(theta), i*b_CC*np.sin(theta), -b_CH*np.cos(theta)])
                graph.add_node(index, name='C')
                graph.add_node(index+1, name='H')
                graph.add_node(index+2, name='H')
                graph.add_edge(index-3, index)
                graph.add_edge(index, index+1)
                graph.add_edge(index, index+2)
                index += 3
                if i == N_carbons-1:
                    pos.append([0, i*b_CC*np.sin(theta)+b_CH*np.sin(theta), b_CH*np.cos(theta)])
                    graph.add_node(index, name='H')
                    graph.add_edge(index-3, index)

        self.N_beads = len(graph)
        self.pos = np.array(pos) * unit.angstrom
        self.graph = graph
        self.b_CH = b_CH
        self.b_CC = b_CC
        assert len(graph) == len(pos), "length of graph != pos"

    def init_topology(self, N_carbons=10):
        top = app.topology.Topology()
        chain = top.addChain('A')
        res = top.addResidue('R1', chain)
        for node in self.graph.nodes:
            element = self.graph.nodes[node]['name']
            top.addAtom(element, app.element.Element._elements_by_symbol[element], res)

        atom_list = list(chain.atoms())

        for e in nx.edge_dfs(self.graph):
            i, j = e[0], e[1]
            top.addBond(atom_list[i], atom_list[j])
        self.top = top

    def init_system(self):
        sys = mm.System()
        sys.addForce(mm.openmm.CMMotionRemover())
        for node in self.graph.nodes:
            element = self.graph.nodes[node]['name']
            sys.addParticle(mass[element])

        self.system = sys

    def write_init_config(self, path='init.pdb'):
        with open(path, 'w') as out:
            app.PDBFile.writeModel(self.top, self.pos, out)
            app.PDBFile.writeFooter(self.top, out)

    #def write_psf(self, path="chain.psf"):
    #    with open(path, 'w') as f:
    #        f.write("{0:8d} !NATOM\n".format(self.top.getNumAtoms()))
    #        for node in self.graph.nodes:
    #            element = self.graph.nodes[node]['name']
    #            f.write('%8d %-4s %-4d %-4s %-4s %-4s %10.6f %13.4f %11d\n' %
    #                    (node+1, "A", 1, "R1", element, element, 0, mass[element], 0))

    #        f.write('\n')
    #        f.write("{0:8d} !NBOND: bonds\n".format(self.top.getNumBonds()))

    #        i = 0
    #        for e in nx.edge_dfs(self.graph):
    #            f.write("%8s%8s" % (e[0], e[1]))
    #            if i % 4 == 0:
    #                f.write("\n")
    #            i += 1

    def add_bond(self):
        bond = mm.openmm.HarmonicBondForce()
        for e in nx.edge_dfs(self.graph):
            element1, element2 = self.graph.nodes[e[0]]['name'], self.graph.nodes[e[1]]['name']
            if element1=='C' and element2=='C':
                bond.addBond(e[0], e[1], self.b_CC*unit.angstrom, 310 * unit.kilocalorie_per_mole / unit.angstrom ** 2)
            else:
                bond.addBond(e[0], e[1], self.b_CH*unit.angstrom, 340 * unit.kilocalorie_per_mole / unit.angstrom ** 2)

        fg_hot = self.name2fg['bond']
        bond.setForceGroup(fg_hot)
        self.fgs.append(fg_hot)
        self.system.addForce(bond)

    def add_angle(self):
        angle = mm.openmm.HarmonicAngleForce()
        for center in self.graph.nodes:
            neighbors = list(self.graph.adj[center])
            for pre, post in combinations(neighbors, 2):
                angle.addAngle(pre, center, post, 109.5 * np.pi / 180 * unit.radians, 50 * unit.kilocalorie_per_mole / unit.radian ** 2)

        fg_hot = self.name2fg['angle']
        angle.setForceGroup(fg_hot)
        self.fgs.append(fg_hot)
        self.system.addForce(angle)

    def add_dihedral(self):
        dihedral = mm.openmm.PeriodicTorsionForce()
        for n1, n2 in nx.edge_dfs(self.graph):
            neighbor1 = list(self.graph.adj[n1])
            neighbor1.remove(n2)
            neighbor2 = list(self.graph.adj[n2])
            neighbor2.remove(n1)
            for pre, post in product(neighbor1, neighbor2):
                if pre != post:
                    dihedral.addTorsion(pre, n1, n2, post, 3, 0 * unit.radian, 0.15 * unit.kilocalorie_per_mole)

        fg_hot = self.name2fg['dihedral']
        dihedral.setForceGroup(fg_hot)
        self.fgs.append(fg_hot)
        self.system.addForce(dihedral)

    def add_LJ(self):
        LJ = mm.openmm.NonbondedForce()
        for node in self.graph.nodes:
            element = self.graph.nodes[node]['name']
            if element == 'C':
                LJ.addParticle(0, 1.9090 * unit.angstrom, 0.1094 * unit.kilocalorie_per_mole)
            else:
                LJ.addParticle(0, 1.4870 * unit.angstrom, 0.0157 * unit.kilocalorie_per_mole)

        bond_indices = [(bond_i.atom1.index, bond_i.atom2.index) for bond_i in self.top.bonds()]
        LJ.createExceptionsFromBonds(bond_indices, 0, 0)

        fg_hot = self.name2fg['LJ']
        LJ.setForceGroup(fg_hot)
        self.fgs.append(fg_hot)
        self.system.addForce(LJ)

    def add_reporter(self, n_steps, report_interval=1000, prefix=".", restart=False):
        self.simulation.reporters.append(app.PDBReporter(f"{prefix}/output.pdb",
                                                         report_interval*10, 
                                                         append=restart,
                                                        ))
        self.simulation.reporters.append(app.DCDReporter(f"{prefix}/output.dcd",
                                                         report_interval,
                                                         append=restart,
                                                        ))
        self.simulation.reporters.append(
                            app.StateDataReporter(f"{prefix}/log", report_interval,
                                                  step=True,
                                                  time=True,
                                                  potentialEnergy=True,
                                                  kineticEnergy=True,
                                                  totalEnergy=True,
                                                  temperature=True,
                                                  elapsedTime=True,
                                                  speed=True,
                                                  remainingTime=True,
                                                  totalSteps=n_steps,
                                                  append=restart,
                                                  ))
        self.simulation.reporters.append(
                        app.CheckpointReporter(f'{prefix}/checkpnt.chk', report_interval))

    def simulate(self, n_steps, n_record, step_size=0.5*unit.femtosecond, restart=False):
        self.system.addForce(mm.AndersenThermostat(298*unit.kelvin, 1/unit.picosecond))
        integrator = mm.VerletIntegrator(step_size)
        self.simulation = app.Simulation(self.top, self.system, integrator)
        self.simulation.context.setPositions(self.pos)
        self.simulation.context.setVelocitiesToTemperature(298*unit.kelvin)
        report_interval = max(1, int(n_steps/n_record))
        self.add_reporter(n_steps, report_interval=report_interval, restart=restart)
        print(f"INFO - We are using {self.simulation.context.getPlatform().getName()} platform.")
        if restart:
            self.simulation.loadCheckpoint('checkpnt.chk')
        self.simulation.step(n_steps)
