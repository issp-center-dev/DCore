import pytriqs.utility.mpi as mpi
from pytriqs.operators.util import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *
from pytriqs.gf.local import *
from pytriqs.applications.dft.sumk_dft import *
from pytriqs.applications.dft.converters.wien2k_converter import *

dft_filename='SrVO3'
U = 9.6
J = 0.8
beta = 40
loops = 10                       # Number of DMFT sc-loops
sigma_mix = 1.0                  # Mixing factor of Sigma after solution of the AIM
delta_mix = 1.0                  # Mixing factor of Delta as input for the AIM
dc_type = 0                      # DC type: 0 FLL, 1 Held, 2 AMF
use_blocks = True                # use bloc structure from DFT input
prec_mu = 0.0001
h_field = 0.0

# Solver parameters
p = {}
p["max_time"] = -1
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 200
p["n_warmup_cycles"] = 100000
p["n_cycles"] = 1000000
p["perfrom_tail_fit"] = True
p["fit_max_moments"] = 4
p["fit_min_n"] = 30
p["fit_max_n"] = 60

# If conversion step was not done, we could do it here. Uncomment the lines it you want to do this.
#from pytriqs.applications.dft.converters.wien2k_converter import *
#Converter = Wien2kConverter(filename=dft_filename, repacking=True)
#Converter.convert_dft_input()
#mpi.barrier()

previous_runs = 0
previous_present = False
if mpi.is_master_node():
    f = HDFArchive(dft_filename+'.h5','a')
    if 'dmft_output' in f:
        ar = f['dmft_output']
        if 'iterations' in ar:
            previous_present = True
            previous_runs = ar['iterations']
    else:
        f.create_group('dmft_output')
    del f
previous_runs    = mpi.bcast(previous_runs)
previous_present = mpi.bcast(previous_present)

SK=SumkDFT(hdf_file=dft_filename+'.h5',use_dft_blocks=use_blocks,h_field=h_field)

n_orb = SK.corr_shells[0]['dim']
l = SK.corr_shells[0]['l']
spin_names = ["up","down"]
orb_names = [i for i in range(n_orb)]

# Use GF structure determined by DFT blocks
gf_struct = SK.gf_struct_solver[0]

# Construct Slater U matrix 
Umat = U_matrix(n_orb=n_orb, U_int=U, J_hund=J, basis='cubic',)

# Construct Hamiltonian and solver
h_int = h_int_slater(spin_names, orb_names, map_operator_structure=SK.sumk_to_solver[0], U_matrix=Umat)
S = Solver(beta=beta, gf_struct=gf_struct)

if previous_present:
  chemical_potential = 0
  dc_imp = 0
  dc_energ = 0
  if mpi.is_master_node():
      ar = HDFArchive(dft_filename+'.h5','a')
      S.Sigma_iw << ar['dmft_output']['Sigma_iw']
      del ar
      chemical_potential,dc_imp,dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])
  S.Sigma_iw << mpi.bcast(S.Sigma_iw)
  chemical_potential = mpi.bcast(chemical_potential)
  dc_imp = mpi.bcast(dc_imp)
  dc_energ = mpi.bcast(dc_energ)
  SK.set_mu(chemical_potential)
  SK.set_dc(dc_imp,dc_energ)

for iteration_number in range(1,loops+1):
    if mpi.is_master_node(): print "Iteration = ", iteration_number

    SK.symm_deg_gf(S.Sigma_iw,orb=0)                        # symmetrise Sigma
    SK.set_Sigma([ S.Sigma_iw ])                            # set Sigma into the SumK class
    chemical_potential = SK.calc_mu( precision = prec_mu )  # find the chemical potential for given density
    S.G_iw << SK.extract_G_loc()[0]                         # calc the local Green function
    mpi.report("Total charge of Gloc : %.6f"%S.G_iw.total_density())

    # Init the DC term and the real part of Sigma, if no previous runs found:
    if (iteration_number==1 and previous_present==False):
        dm = S.G_iw.density()
        SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)
        S.Sigma_iw << SK.dc_imp[0]['up'][0,0]

    # Calculate new G0_iw to input into the solver:
    if mpi.is_master_node():
        # We can do a mixing of Delta in order to stabilize the DMFT iterations:
        S.G0_iw << S.Sigma_iw + inverse(S.G_iw)
        ar = HDFArchive(dft_filename+'.h5','a')
        if (iteration_number>1 or previous_present):
            mpi.report("Mixing input Delta with factor %s"%delta_mix)
            Delta = (delta_mix * delta(S.G0_iw)) + (1.0-delta_mix) * ar['dmft_output']['Delta_iw']
            S.G0_iw << S.G0_iw + delta(S.G0_iw) - Delta
        ar['dmft_output']['Delta_iw'] = delta(S.G0_iw)
        S.G0_iw << inverse(S.G0_iw)
        del ar

    S.G0_iw << mpi.bcast(S.G0_iw)

    # Solve the impurity problem:
    S.solve(h_int=h_int, **p)

    # Solved. Now do post-processing:
    mpi.report("Total charge of impurity problem : %.6f"%S.G_iw.total_density())

    # Now mix Sigma and G with factor sigma_mix, if wanted:
    if (iteration_number>1 or previous_present):
        if mpi.is_master_node():
            ar = HDFArchive(dft_filename+'.h5','a')
            mpi.report("Mixing Sigma and G with factor %s"%sigma_mix)
            S.Sigma_iw << sigma_mix * S.Sigma_iw + (1.0-sigma_mix) * ar['dmft_output']['Sigma_iw']
            S.G_iw << sigma_mix * S.G_iw + (1.0-sigma_mix) * ar['dmft_output']['G_iw']
            del ar
        S.G_iw << mpi.bcast(S.G_iw)
        S.Sigma_iw << mpi.bcast(S.Sigma_iw)

    # Write the final Sigma and G to the hdf5 archive:
    if mpi.is_master_node():
        ar = HDFArchive(dft_filename+'.h5','a')
        ar['dmft_output']['iterations'] = iteration_number + previous_runs
        ar['dmft_output']['G_tau'] = S.G_tau
        ar['dmft_output']['G_iw'] = S.G_iw
        ar['dmft_output']['Sigma_iw'] = S.Sigma_iw
        ar['dmft_output']['G0-%s'%(iteration_number)] = S.G0_iw
        ar['dmft_output']['G-%s'%(iteration_number)] = S.G_iw
        ar['dmft_output']['Sigma-%s'%(iteration_number)] = S.Sigma_iw
        del ar

    # Set the new double counting:
    dm = S.G_iw.density() # compute the density matrix of the impurity problem
    SK.calc_dc(dm, U_interact = U, J_hund = J, orb = 0, use_dc_formula = dc_type)

    # Save stuff into the dft_output group of hdf5 archive in case of rerun:
    SK.save(['chemical_potential','dc_imp','dc_energ'])
