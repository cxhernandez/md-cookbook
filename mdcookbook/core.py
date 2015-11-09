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
except:
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


def rmNV(s):
    s.strip([atom.name == 'NV' for atom in s.atoms])
    return s


def unpack_pose(pose, noNV=True):
    if load_rosetta:
        s = load_rosetta(pose)
        if noNV:
            s = rmNV(s)
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


def solvate(positions, topology, forcefield, ion_content, boxSize=None,
            padding=None, model='tip3p'):
    modeller = Modeller(topology, positions)
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        padding=padding *
        nanometers,
        boxSize=boxSize,
        positiveIon='Na+',
        negativeIon='Cl-',
        ionicStrength=ion_content *
        molar)
    return modeller


def smart_solvate(positions, topology, forcefield, ion_content, n,
                  model='tip3p', tries=10):
    """ Solvate a system with a fixed number of waters

        Parameters
        ----------
        positions : OpenMM :attr:`positions`
        topology : OpenMM :attr:`topology`
        forcefield : OpenMM :class:`ForceField`
        ion_content : Molar concentration of solvent ions (float)
        n : Target number of waters to solvate system (int)
        model : Water model (Default: 'tip3p')
        tries : Number of attempts to reach taget number of waters
                (Default: 10)

        Returns
        -------
        modeller : OpenMM :class:`Modeller`
    """

    # Get initial estimates of the box volume
    modeller = solvate(positions, topology, forcefield, ion_content,
                       model=model, padding=0.0)
    box_o = modeller.topology.getUnitCellDimensions()
    n_wat_o = get_num_res(modeller.topology)
    volume_o = box_o[0] * box_o[1] * box_o[2]

    if n_wat_o > n:
        raise Exception("Target number of waters is too small.")

    # Slowly increase the box size until just above target number of waters
    scale = 0.9 * (n / n_wat_o) ** (1.0 / 3.0)
    over_target = False
    xwat = int(.01 * n_wat_o)
    density = None
    while not over_target and tries > 0:
        modeller = solvate(positions, topology, forcefield, ion_content,
                           model=model, boxSize=scale * box_o)
        n_wat = get_num_res(modeller.topology)
        if (n_wat > n):
            over_target = True
        else:
            if density is None:
                box = modeller.topology.getUnitCellDimensions()
                volume = box[0] * box[1] * box[2]
                density = (n_wat - n_wat_o) / (volume - volume_o)
            delta = (n + xwat - n_wat_o) / density
            scale = ((volume_o + delta) / volume_o) ** (1.0 / 3.0)
            xwat += xwat
            tries -= 1

    # Delete waters to achieve target number
    n_wat_del = n_wat - n
    if n_wat_del > 0:
        modeller = del_res(modeller, n_wat_del)

    if get_num_res(modeller.topology) != n:
        raise Exception("Target solvation could not be completed "
                        "in %d tries." % tries)

    # Adjust ion concentrations to expected concentration
    n_anion = get_num_res(modeller.topology, res_type='CL')
    n_anion_del = int(float(n_wat_del) / n_wat * n_anion)
    n_cation = get_num_res(modeller.topology, res_type='NA')
    n_cation_del = int(float(n_wat_del) / n_wat * n_cation)
    if n_anion_del > 0:
        modeller = del_res(modeller, n_wat_del, res_type='CL')
    if n_cation_del > 0:
        modeller = del_res(modeller, n_wat_del, res_type='NA')

    return modeller
