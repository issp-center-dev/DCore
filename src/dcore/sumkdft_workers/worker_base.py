from dcore._dispatcher import mpi
#from dcore._dispatcher import GfReFreq, GfImFreq
from dcore._dispatcher import HDFArchive

# Run SumkDFT using MPI

def _add_potential(_sigma, _pot):
    sigma_plus_pot = _sigma.copy()
    for sp, sigma in sigma_plus_pot:
        sigma += _pot[sp]
    return sigma_plus_pot
    

def setup_sk(sk, iwn_or_w_or_none, params):
    if iwn_or_w_or_none == 'iwn':
        assert len(params['Sigma_iw_sh']) == len(params['potential'])
        Sigma_iw_sh_plus_pot = [_add_potential(sigma, pot)
                                for sigma, pot in zip(params['Sigma_iw_sh'], params['potential'])]
        sk.set_Sigma(Sigma_iw_sh_plus_pot)
    elif iwn_or_w_or_none == 'w':
        # sk.set_Sigma([params['Sigma_w_sh'][ish] for ish in range(sk.n_inequiv_shells)])
        Sigma_w_sh = [params['Sigma_w_sh'][ish] for ish in range(sk.n_inequiv_shells)]
        Sigma_w_sh_plus_pot = [_add_potential(sigma, pot)
                           for sigma, pot in zip(Sigma_w_sh, params['potential'])]
        sk.set_Sigma(Sigma_w_sh_plus_pot)
    elif iwn_or_w_or_none == "none":
        pass
    else:
        raise RuntimeError("Invalid iwn_or_w")
    
    if params['with_dc']:
        sk.set_dc(params['dc_imp'], params['dc_energ'])
    sk.set_mu(params['mu'])


class SumkDFTWorkerBase(object):
    """
    Perform some MPI calculations with a subclass of SumkDFT
    """
    def __init__(self, model_hdf5_file, input_file, output_file) -> None:
        """
        model_hdft_file: str
            Model file name
        input_file: str
            Name of HDF5 file that contains input data
        output_file: str
            Name of HDF5 file to which results are stored.
        """
        super().__init__()
        self.model_hdf5_file = model_hdf5_file
        self.input_file = input_file
        self.output_file = output_file

        # read HDF5 file on the master node
        if mpi.is_master_node():
            with HDFArchive(input_file, 'r') as h:
                params = h['params']
                keys = list(params.keys())
        else:
            params = {}
            keys = []
        assert isinstance(params, dict)
        self.params = params
    
        # broadcast parameters
        keys = mpi.bcast(keys)
        for key in keys:
            if not mpi.is_master_node():
                params[key] = None
            params[key] = mpi.bcast(params[key])
    
    
    def run(self):
        raise NotImplementedError()

    def save_result(self, results: dict):
        if not mpi.is_master_node():
            raise RuntimeError("save_result must be called on the master node!")

        with HDFArchive(self.output_file, 'w') as h:
            for k, v in list(results.items()):
                h[k] = v