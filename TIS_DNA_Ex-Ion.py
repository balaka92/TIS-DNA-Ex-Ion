import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import numpy as np
import math
import os, sys
import json
import datetime
from sys import platform
from sys import stdout


#********************************OPENMM CODE TO SIMULATE SSDNA********************************#


###################~~~~~~~~Keep a track of simulation time~~~~~~~~###################

time_st = datetime.datetime.utcnow()

###################~~~~~~~~Load the input file~~~~~~~~###################

configdic = json.load(open(sys.argv[1]))

###################~~~~~~~~Print the OPENMM version~~~~~~~~###################

version=mm.Platform.getOpenMMVersion()
print('simulations are running on OPENMM\t')
print(version)

###################~~~~~~~~Set input and output file names~~~~~~~~###################

datafile = open('output.dat', 'a')
###################~~~~~~~~Check for restart~~~~~~~~###################

if 'rstin_prefix' in configdic.keys():
    rstin_prefix = configdic['rstin_prefix']
if 'rstout_prefix' in configdic.keys():
    rstout_prefix = configdic['rstout_prefix']

###################~~~~~~~~Read system composition from input file~~~~~~~~###################

Nnuc=configdic['Nnuc']
temp=configdic['temp']
Box_len=configdic['Box_len']



###################~~~~~~~~Set simulation length and data interval frequency~~~~~~~~###################

numsteps=configdic['numsteps']
data_interval=configdic['data_interval']

###################~~~~~~~~Set up the simulation~~~~~~~~###################

pdb = app.PDBFile('INIT.pdb')
positions = pdb.positions
system = mm.System()
system.setDefaultPeriodicBoxVectors(mm.Vec3(Box_len*u.angstroms, 0, 0), mm.Vec3(0, Box_len*u.angstroms, 0), mm.Vec3(0, 0, Box_len*u.angstroms))

with open('Mass.dat') as input_file:
    for line in input_file:
       columns = line.split()
       mass_atom = float(columns[0])
       system.addParticle(mass_atom*u.amu)


###################~~~~~~~~Read the list of bonds and add the potential~~~~~~~~###################


force_bonded = mm.HarmonicBondForce()
force_bonded.setUsesPeriodicBoundaryConditions(True)
with open('BONDED.dat') as input_file:
    for line in input_file:
        columns = line.split()
        bond_atom_index_i = int(columns[0])
        bond_atom_index_j = int(columns[1])
        bond_k = 2.0*float(columns[2])
        bond_r0 = float(columns[3])
        force_bonded.addBond(bond_atom_index_i, bond_atom_index_j, bond_r0*u.angstroms, bond_k*u.kilocalories_per_mole/u.angstroms**2)
system.addForce(force_bonded)


###################~~~~~~~~Read the list of angles and add the potential~~~~~~~~###################

force_angle = mm.HarmonicAngleForce()
force_angle.setUsesPeriodicBoundaryConditions(True)
with open('ANG.dat') as input_file:
   for line in input_file:
       angle_columns = line.split()
       ang_atom_index_i = int(angle_columns[0])
       ang_atom_index_j = int(angle_columns[1]) 
       ang_atom_index_k = int(angle_columns[2])
       ang_k = 2.0*float(angle_columns[3])
       ang_theta = float(angle_columns[4])
       force_angle.addAngle(ang_atom_index_i, ang_atom_index_j, ang_atom_index_k, ang_theta*u.radians, ang_k*u.kilocalories_per_mole/u.radians**2)
system.addForce(force_angle)


 
###################~~~~~~~~Add Electrostatics using PME~~~~~~~~###################

##READ1: Particle Charges are rescaled by sqrt(water dielectric constant): q_scaled = q /sqrt(water dielectric constant)##
##Default LJ Interactions are set to 0 ##


temp_cent = temp - 273.15
scale_factor = 87.74-(0.4008*temp_cent)+(0.0009398*temp_cent*temp_cent)-(1.41*temp_cent*temp_cent*temp_cent/1000000) ;
CoulombOpenmm = mm.NonbondedForce()
CoulombOpenmm.setNonbondedMethod(mm.NonbondedForce.PME)
CoulombOpenmm.setCutoffDistance(0.3*Box_len*u.angstroms)
CoulombOpenmm.setEwaldErrorTolerance(0.005)
#CoulombOpenmm.setUseDispersionCorrection(True)
with open('Coulomb.dat') as input_file:
   for line in input_file:
      charge_columns = line.split()
      q = float(charge_columns[0])/math.sqrt(scale_factor)
      sig = 0.0
      eps = 0.0
      CoulombOpenmm.addParticle(q*u.elementary_charge, sig*u.angstroms, eps*u.kilocalories_per_mole)
with open('Exclusion_list.dat') as input_file:
   for line in input_file:
      exclu_columns = line.split()
      index1 = int(exclu_columns[0])
      index2 = int(exclu_columns[1])
      CoulombOpenmm.addException(index1, index2, 0.0, 0.0, 0.0)
system.addForce(CoulombOpenmm)
MTHD=CoulombOpenmm.getNonbondedMethod()
CHK=CoulombOpenmm.usesPeriodicBoundaryConditions()




###################~~~~~~~~Add excluded volume interactions~~~~~~~~###################

##############~~~Between DNA and ions~~~###########

ExVol1Force = mm.CustomNonbondedForce( '(scale_factor/2.0)*eps*((val/(val + 0.5*(r-sig)*(1 - ((r-sig)/abs(r-sig)))))^12 - 2.0*(val/(val+0.5*(r-sig)*(1 - ((r-sig)/abs(r-sig)))))^6 + 1.0);  scale_factor = abs(scale_factor1 - scale_factor2); eps=sqrt(eps1*eps2); sig=(sig1 + sig2)')
ExVol1Force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
ExVol1Force.addPerParticleParameter('scale_factor')
ExVol1Force.addPerParticleParameter('eps')
ExVol1Force.addPerParticleParameter('sig')
ExVol1Force.addGlobalParameter("val", 1.6*u.angstroms)
ExVol1Force.setCutoffDistance(6.0*u.angstroms)

with open('EV_DNA_ION.dat') as input_file:
   for line in input_file:
      ev1_columns = line.split()
      sig = float(ev1_columns[0])
      eps = float(ev1_columns[1])
      scale_factor = float(ev1_columns[2])
      ExVol1Force.addParticle([scale_factor, eps*u.kilocalories_per_mole, sig*u.angstroms])
with open('Exclusion_list.dat') as input_file:
   for line in input_file:
      exclu_columns = line.split()
      index1 = int(exclu_columns[0])
      index2 = int(exclu_columns[1])
      ExVol1Force.addExclusion(index1, index2)
system.addForce(ExVol1Force)

#############~~~Between ions and ions~~~##################

ExVol2Force = mm.CustomNonbondedForce( 'eps*((val/(val + 0.5*(r-sig)*(1 - ((r-sig)/abs(r-sig)))))^12 - 2.0*(val/(val+0.5*(r-sig)*(1 - ((r-sig)/abs(r-sig)))))^6 + 1.0);  eps=sqrt(eps1*eps2); sig=(sig1 + sig2)')
ExVol2Force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
ExVol2Force.addPerParticleParameter('eps')
ExVol2Force.addPerParticleParameter('sig')
ExVol2Force.addGlobalParameter("val", 1.6*u.angstroms)
ExVol2Force.setCutoffDistance(6.0*u.angstroms)

with open('EV_ION_ION.dat') as input_file:
   for line in input_file:
      ev2_columns = line.split()
      sig = float(ev2_columns[0])
      eps = float(ev2_columns[1])
      ExVol2Force.addParticle([eps*u.kilocalories_per_mole, sig*u.angstroms])
with open('Exclusion_list.dat') as input_file:
   for line in input_file:
      exclu_columns = line.split()
      index1 = int(exclu_columns[0])
      index2 = int(exclu_columns[1])
      ExVol2Force.addExclusion(index1, index2)
system.addForce(ExVol2Force)

###############~~~Between DNA beads~~~##################

ExVol3Force = mm.CustomNonbondedForce( 'eps*((val/(r - val))^12 - 2.0*(val/(r - val))^6 + 1.0);  eps=sqrt(eps1*eps2);')
ExVol3Force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
ExVol3Force.addPerParticleParameter('eps')
ExVol3Force.addGlobalParameter("val", 1.6*u.angstroms)
ExVol3Force.setCutoffDistance(3.2*u.angstroms)

with open('EV_DNA_DNA.dat') as input_file:
   for line in input_file:
      ev3_columns = line.split()
      sig = float(ev3_columns[0])
      eps = float(ev3_columns[1])
      ExVol3Force.addParticle([eps*u.kilocalories_per_mole])
with open('Exclusion_list.dat') as input_file:
   for line in input_file:
      exclu_columns = line.split()
      index1 = int(exclu_columns[0])
      index2 = int(exclu_columns[1])
      ExVol3Force.addExclusion(index1, index2)
system.addForce(ExVol3Force) 


###################~~~~~~~~Add stacking interactions~~~~~~~~###################     


Stackingforce = mm.CustomCompoundBondForce(7, "U0/(1.0 + kbond*(distance(p6,p7) - r0)^2  + kphi1*(dihedral(p1,p2,p3,p4) - phi10 + pi*(((pi - dihedral(p1,p2,p3,p4) + phi10)/abs(pi - dihedral(p1,p2,p3,p4) + phi10)) -  ((pi + dihedral(p1,p2,p3,p4) - phi10)/abs(pi + dihedral(p1,p2,p3,p4) - phi10)))  )^2  +  kphi2*(dihedral(p2,p3,p4,p5) - phi20 + pi*(((pi - dihedral(p2,p3,p4,p5) + phi20)/abs(pi - dihedral(p2,p3,p4,p5) + phi20)) -  ((pi + dihedral(p2,p3,p4,p5) - phi20)/abs(pi + dihedral(p2,p3,p4,p5) - phi20)))  )^2 )");
Stackingforce.addPerBondParameter("U0");
Stackingforce.addPerBondParameter("r0");
Stackingforce.addPerBondParameter("phi10");
Stackingforce.addPerBondParameter("phi20");
Stackingforce.setUsesPeriodicBoundaryConditions(True)
Stackingforce.addGlobalParameter('kbond', 1.45/u.angstroms**2)
Stackingforce.addGlobalParameter('kphi1', 3.0/u.radians**2)
Stackingforce.addGlobalParameter('kphi2', 3.0/u.radians**2)
Stackingforce.addGlobalParameter('pi', np.pi)
Stackingforce.setUsesPeriodicBoundaryConditions(True)

with open('STACK.dat') as input_file:
   for line in input_file:
      stack_columns = line.split()
      N1 = int(stack_columns[0])
      N2 = int(stack_columns[1])
      r0 = float(stack_columns[2])
      phi10 = float(stack_columns[3])
      phi20 = float(stack_columns[4])
      U0 = float(stack_columns[5])
      p1=3*N1-1+N2;
      p2=3*N1+N2;
      p3=3*N1+2+N2;
      p4=3*N1+3+N2;
      p5=3*N1+5+N2;
      p6=3*N1+1+N2;
      p7=3*N1 + 4+N2;
      group_add = [p1, p2, p3, p4, p5, p6, p7]
      Stackingforce.addBond(group_add, [U0*u.kilocalories_per_mole,r0*u.angstroms, phi10*u.radians, phi20*u.radians])
system.addForce(Stackingforce)


###################~~~~~~~~Add hydrogen bonding interactions~~~~~~~~###################


Hbondingforce = mm.CustomCompoundBondForce(6,"scale1*UHyd/(1.0 + kdist*DIST*DIST + kang*ANG1*ANG1 + kang*ANG2*ANG2 + kdih*DIHED1*DIHED1 + kdih*DIHED2*DIHED2 + kdih*DIHED3*DIHED3); scale1 = 0.25*(1.0 - ((angle(p2,p1,p4)-pi_cut)/(abs(angle(p2,p1,p4)-pi_cut))))*(1.0 - ((angle(p5,p4,p1)-pi_cut)/(abs(angle(p5,p4,p1)-pi_cut)))); DIST = (distance(p1,p4) - dist_0); ANG1 = (angle(p2,p1,p4) - angle_01); ANG2 = (angle(p5,p4,p1) - angle_02); DIHED1 = dihedral(p2,p1,p4,p5) - dihedral_01 + pi*(((pi - dihedral(p2,p1,p4,p5) + dihedral_01)/abs(pi - dihedral(p2,p1,p4,p5) + dihedral_01)) - ((pi + dihedral(p2,p1,p4,p5) - dihedral_01)/abs(pi + dihedral(p2,p1,p4,p5) - dihedral_01))) ;  DIHED2 = (dihedral(p6,p5,p4,p1) - dihedral_02 + pi*(((pi - dihedral(p6,p5,p4,p1) + dihedral_02)/abs(pi - dihedral(p6,p5,p4,p1) + dihedral_02)) - ((pi + dihedral(p6,p5,p4,p1) - dihedral_02)/abs(pi + dihedral(p6,p5,p4,p1) - dihedral_02)))); DIHED3 = (dihedral(p3,p2,p1,p4) - dihedral_03 + pi*(((pi - dihedral(p3,p2,p1,p4) + dihedral_03)/abs(pi - dihedral(p3,p2,p1,p4) + dihedral_03)) - ((pi + dihedral(p3,p2,p1,p4) - dihedral_03)/abs(pi + dihedral(p3,p2,p1,p4) - dihedral_03))));")


Hbondingforce.addPerBondParameter("UHyd")
Hbondingforce.addPerBondParameter("dist_0")
Hbondingforce.addPerBondParameter("angle_01")
Hbondingforce.addPerBondParameter("angle_02")
Hbondingforce.addPerBondParameter("dihedral_01")
Hbondingforce.addPerBondParameter("dihedral_02")
Hbondingforce.addPerBondParameter("dihedral_03")
Hbondingforce.addGlobalParameter('pi_cut',3.124139361*u.radians)
Hbondingforce.addGlobalParameter('kdist', 4.00/u.angstroms**2)
Hbondingforce.addGlobalParameter('kang', 1.50/u.radians**2)
Hbondingforce.addGlobalParameter('kdih', 0.15/u.radians**2)
Hbondingforce.addGlobalParameter('pi', np.pi)
Hbondingforce.setUsesPeriodicBoundaryConditions(True)

with open('HYD_BOND.dat') as input_file:
   for line in input_file:
      hbond_columns = line.split()
      p1 = int(hbond_columns[0])
      p2 = int(hbond_columns[1])
      p3 = int(hbond_columns[2])
      p4 = int(hbond_columns[3])
      p5 = int(hbond_columns[4])
      p6 = int(hbond_columns[5])
      dist_0 = float(hbond_columns[6])
      angle_01 = float(hbond_columns[7])
      angle_02 = float(hbond_columns[8])
      dihedral_01 = float(hbond_columns[9])
      dihedral_02 = float(hbond_columns[10])
      dihedral_03 = float(hbond_columns[11])
      UHyd = float(hbond_columns[12])
      group_add2 = [p1,p2,p3,p4,p5,p6]
      Hbondingforce.addBond(group_add2, [UHyd*u.kilocalories_per_mole, dist_0*u.angstroms, angle_01*u.radians, angle_02*u.radians, dihedral_01*u.radians, dihedral_02*u.radians, dihedral_03*u.radians])
system.addForce(Hbondingforce)

###################~~~~~~~~Compare energy to check with my simulation code~~~~~~~~###################


############Group the forces###############

for i in range(system.getNumForces()):
    force = system.getForce(i)
    force.setForceGroup(i)


###################~~~~~~~~Pick the integrator~~~~~~~~###################

integrator = mm.LangevinIntegrator(temp*u.kelvin, 0.01/u.picosecond, 2.0*u.femtosecond)


###################~~~~~~~~Simulation block~~~~~~~~###################

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'Precision': 'double'}
simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temp*u.kelvin)
simulation.minimizeEnergy(maxIterations=25000)

###################~~~~~~~~If restart is applicable then read the position and velocities~~~~~~~~###################
restart=False
if 'rstin_prefix' in configdic.keys():
    restart=True
    with open('checkpnt.chk', 'rb') as f:
        simulation.context.loadCheckpoint(f.read())


###################~~~~~~~~Print the energy decomposition~~~~~~~~###################

print(simulation.context.getState(getEnergy=True,groups={0}).getPotentialEnergy())
print(simulation.context.getState(getEnergy=True,groups={1}).getPotentialEnergy())
print(simulation.context.getState(getEnergy=True,groups={2}).getPotentialEnergy())
print(simulation.context.getState(getEnergy=True,groups={3}).getPotentialEnergy())
print(simulation.context.getState(getEnergy=True,groups={4}).getPotentialEnergy())
print(simulation.context.getState(getEnergy=True,groups={5}).getPotentialEnergy())
print(simulation.context.getState(getEnergy=True,groups={6}).getPotentialEnergy())
print(simulation.context.getState(getEnergy=True,groups={7}).getPotentialEnergy())


###################~~~~~~~~Report the simulation progress~~~~~~~~###################

simulation.reporters.append(app.DCDReporter('ds.dcd', data_interval, append = restart)) 
simulation.reporters.append(app.StateDataReporter(stdout, data_interval, step=True, potentialEnergy=True,  temperature=True))
simulation.reporters.append(app.StateDataReporter(datafile, data_interval, step=True, potentialEnergy=True,  temperature=True))
simulation.reporters.append(app.CheckpointReporter('checkpnt.chk', data_interval))
simulation.step(numsteps)



###################~~~~~~~~Report the time taken by the simulation~~~~~~~~###################

time_en = datetime.datetime.utcnow()
time_lapse = time_en - time_st
print('Simulation took\t')
print(time_lapse)
