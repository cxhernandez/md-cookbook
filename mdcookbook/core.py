from simtk.openmm import LangevinIntegrator, MonteCarloBarostat
from simtk.openmm.app import ForceField, Simulation, Modeller, HBonds, PME
from simtk.unit import picoseconds, femtoseconds, nanometers
from simtk.unit import kelvin, molar, bar

from random import shuffle
from mdcookbook.utils import count

try:
    from rosetta import Pose, FastRelax
    from rosetta import hbond_lr_bb, fa_pair, fa_elec, ref, rama
    from rosetta import pose_from_pdb, pose_from_sequence, get_fa_scorefxn
    from toolbox import mutate_residue as mutate

    TERMS = [hbond_lr_bb, fa_pair, fa_elec, ref, rama]

except ImportError:
    FastRelax = get_fa_scorefxn = pose_from_pdb = mutate = None

try:
    from chemistry import load_rosetta
except ImportError:
    load_rosetta = None


def add_caps(pose):
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


def model_from_pdb(pdb, mut_pos=None, mut=None, cap=True):
    if pose_from_pdb is None:
        raise ImportError('Could not find PyRosetta.')
    pose = pose_from_pdb(pdb)
    scorefxn = get_fa_scorefxn()
    relax = FastRelax(scorefxn, 15)
    if cap:
        pose = add_caps(pose)
    if mut_pos:
        pose = mutate(pose, mut_pos, mut)
    relax.apply(pose)
    return pose


def model_from_seq(seq, cap=True):
    if pose_from_sequence is None:
        raise ImportError('Could not find PyRosetta.')
    pose = pose_from_sequence(seq)
    scorefxn = get_fa_scorefxn()
    for term in TERMS:
        scorefxn.set_weight(term, 1.0)
    relax = FastRelax(scorefxn, 30)
    if cap:
        pose = add_caps(pose)
    relax.apply(pose)
    return pose


def get_sim(positions, topology, temp, forcefield, platform, props, nbcutoff=1,
            intstep=2, baro=1, tol=1e-5):
    system = forcefield.createSystem(topology, nonbondedMethod=PME,
                                     nonbondedCutoff=nbcutoff * nanometers,
                                     constraints=HBonds)
    integrator = LangevinIntegrator(temp * kelvin, 1 / picoseconds,
                                    intstep * femtoseconds)
    integrator.setConstraintTolerance(tol)
    system.addForce(MonteCarloBarostat(baro * bar, temp * kelvin))
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


def rmnv(s):
    s.strip([atom.name == 'NV' for atom in s.atoms])
    return s


def unpack_pose(pose, nonv=True):
    if load_rosetta:
        s = load_rosetta(pose)
        if nonv:
            s = rmnv(s)
    else:
        raise ImportError('Could not find ParmEd.')
    return s.positions, s.topology


def get_res(topology, res_type='HOH'):
    for residue in topology.residues():
        if residue.name == res_type:
            yield residue


def get_num_res(topology, res_type='HOH'):
    return count(get_res(topology, res_type))


def del_res(modeller, n_del, res_type='HOH'):
    res = list(get_res(modeller.topology, res_type=res_type))
    shuffle(res)
    modeller.delete(res[:n_del])
    return modeller


def solvate(positions, topology, forcefield, ion_content, model='tip3p',
            neutralize=True, **kwargs):
    modeller = Modeller(topology, positions)
    modeller.addSolvent(
        forcefield,
        model=model,
        positiveIon='Na+',
        negativeIon='Cl-',
        ionicStrength=ion_content*molar,
        neutralize=neutralize,
        **kwargs)
    return modeller
