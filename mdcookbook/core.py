from simtk.openmm import LangevinIntegrator, MonteCarloBarostat
from simtk.openmm.app import ForceField, Simulation, Modeller, HBonds, PME
from simtk.unit import picoseconds, femtoseconds, nanometers
from simtk.unit import kelvin, molar, bar

try:
    from rosetta import Pose, pose_from_pdb, get_fa_scorefxn, FastRelax
    from toolbox import mutate_residue as mutate
except ImportError:
    FastRelax = get_fa_scorefxn = pose_from_pdb = mutate = None

try:
    from chemistry import load_rosetta
except:
    load_rosetta = None


def addCaps(pose):
    if mutate is None:
        raise ImportError('Could not find PyRosetta.')
    if not isinstance(pose, Pose):
        raise TypeError("pose should be a PyRosetta Pose object.")
    first_rsym = pose.sequence()[0]
    first_rname = pose.residue(1).name3().strip()
    last_rsym = pose.sequence()[-1]
    last_rname = pose.residue(pose.n_residue()).name3().strip()

    pose = mutate(pose, 1,
                  '%s[%s_p:N_acetylated]' % (first_rsym, first_rname))
    return mutate(pose, pose.n_residue(),
                  '%s[%s_p:C_methylamidated]' % (last_rsym, last_rname))


def model(pdb, mutpos=None, mut=None, cap=True):
    if mutate is None:
        raise ImportError('Could not find PyRosetta.')
    scorefxn = get_fa_scorefxn()
    relax = FastRelax(scorefxn, 5)
    pose = pose_from_pdb(pdb)
    if cap:
        pose = addCaps(pose)
    if mutpos:
        pose = mutate(pose, mutpos, mut)
    relax.apply(pose)
    return pose


def solvate(positions, topology, ff, smolar, boxsize):
    mod = Modeller(topology, positions)
    mod.addSolvent(ff, model='tip3p', boxSize=boxsize*nanometers,
                   positiveIon='Na+', negativeIon='Cl-',
                   ionicStrength=smolar*molar)
    return mod


def get_sim(positions, topology, temp, ff, platform, props, nbcutoff=1,
            intstep=2, baro=1, tol=1e-5):
    system = ff.createSystem(topology, nonbondedMethod=PME,
                             nonbondedCutoff=nbcutoff*nanometers,
                             constraints=HBonds)
    integrator = LangevinIntegrator(temp*kelvin, 1/picoseconds,
                                    intstep*femtoseconds)
    integrator.setConstraintTolerance(tol)
    system.addForce(MonteCarloBarostat(baro*bar, temp*kelvin))
    simulation = Simulation(topology, system, integrator, platform, props)
    simulation.context.setPositions(positions)
    return simulation, system, integrator


def get_state(simulation, temp):
    simulation.context.setVelocitiesToTemperature(temp)
    return simulation.context.getState(getPositions=True, getVelocities=True,
                                       getForces=True, getEnergy=True,
                                       getParameters=True,
                                       enforcePeriodicBox=True)


def get_ff(param='amber99sbildn.xml', watmod='tip3p.xml'):
    return ForceField(param, watmod)


def unpack_pose(pose):
    if load_rosetta:
        s = load_rosetta(pose)
    else:
        raise ImportError('Could not find ParmEd.')
    return s.positions, s.topology
