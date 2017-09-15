from pytriqs.applications.dft.sumk_dft_tools import *
from pytriqs.applications.dft.converters.wien2k_converter import *
from pytriqs.applications.impurity_solvers.hubbard_I.hubbard_solver import Solver

# Creates the data directory, cd into it:
#Prepare_Run_Directory(DirectoryName = "Ce-Gamma") 
dft_filename = 'Ce-gamma'
beta =  40
U_int = 6.00
J_hund = 0.70
DC_type = 0                      # 0...FLL, 1...Held, 2... AMF, 3...Lichtenstein
ommin=-4.0
ommax=6.0
N_om=2001
broadening = 0.02

# Convert DMFT input:
Converter = Wien2kConverter(filename=dft_filename,repacking=True)
Converter.convert_dft_input()
Converter.convert_parproj_input()

# Init the SumK class
SK = SumkDFTTools(hdf_file=dft_filename+'.h5',use_dft_blocks=False)

# load old chemical potential and DC
if mpi.is_master_node():
    SK.chemical_potential,SK.dc_imp,SK.dc_energ = SK.load(['chemical_potential','dc_imp','dc_energ'])

SK.chemical_potential = mpi.bcast(SK.chemical_potential)
SK.dc_imp = mpi.bcast(SK.dc_imp)
SK.dc_energ = mpi.bcast(SK.dc_energ)

if (mpi.is_master_node()):
    print 'DC after reading SK: ',SK.dc_imp[0]

N = SK.corr_shells[0]['dim']
l = SK.corr_shells[0]['l']

# Init the Solver:
S = Solver(beta = beta, l = l)

# set atomic levels:
eal = SK.eff_atomic_levels()[0]
S.set_atomic_levels( eal = eal )

# Run the solver to get GF and self-energy on the real axis
S.GF_realomega(ommin=ommin, ommax = ommax, N_om=N_om,U_int=U_int,J_hund=J_hund)
SK.set_Sigma([S.Sigma])

# compute DOS
SK.dos_parproj_basis(broadening=broadening)
