from dcore._dispatcher import mpi
from .worker_base import SumkDFTWorkerBase, setup_sk

class SumkDFTWorkerGloc(SumkDFTWorkerBase):
    """For computing Gloc"""
    def __init__(self, model_hdf5_file, input_file, output_file) -> None:
        super().__init__(model_hdf5_file, input_file, output_file)

    def run(self):
        from dcore.sumkdft_opt import SumkDFT_opt as SumkDFT
        beta = self.params['beta']
        with_dc = self.params['with_dc']
        results = {}

        sk = SumkDFT(hdf_file=self.model_hdf5_file, use_dft_blocks=False, h_field=0.0)
        setup_sk(sk, 'iwn', self.params)

        # workaround for 'High frequency error'
        if self.params['no_tail_fit']:
            sk.total_density = sk.total_density_matsubara
            sk.density_matrix = sk.density_matrix_matsubara

        if self.params['adjust_mu']:
            # find the chemical potential for given density
            sk.calc_mu(self.params['prec_mu'])
            if mpi.is_master_node():
                results['mu'] = float(sk.chemical_potential)

        # Local Green's function and Density matrix
        Gloc = sk.extract_G_loc(with_dc=with_dc)
        dm_corr_sh = sk.density_matrix(beta=beta)
        dm_sh = [dm_corr_sh[sk.inequiv_to_corr[ish]] for ish in range(sk.n_inequiv_shells)]
        if mpi.is_master_node():
            results['Gloc_iw_sh'] = Gloc
            results['dm_sh'] = dm_sh
            self.save_result(results)
