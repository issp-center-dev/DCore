import triqs.utility.mpi as mpi
from .worker_base import SumkDFTWorkerBase, setup_sk

class SumkDFTWorkerMomDist(SumkDFTWorkerBase):
    """For computing momentum_distribution"""
    def __init__(self, model_hdf5_file, input_file, output_file) -> None:
        super().__init__(model_hdf5_file, input_file, output_file)

    def run(self):
        from dcore.sumkdft_post import SumkDFTDCorePost
        sk = SumkDFTDCorePost(hdf_file=self.model_hdf5_file, use_dft_blocks=False, h_field=0.0)
        setup_sk(sk, 'iwn', self.params)
        results = {}
        results['den'] = \
            sk.calc_momentum_distribution(mu=self.params["mu"], beta=self.params['beta'],
            with_Sigma=True, with_dc=True)

        if mpi.is_master_node():
            self.save_result(results)