import triqs.utility.mpi as mpi
from .worker_base import SumkDFTWorkerBase, setup_sk

class SumkDFTWorkerDOS(SumkDFTWorkerBase):
    """For computing DOS"""
    def __init__(self, model_hdf5_file, input_file, output_file) -> None:
        super().__init__(model_hdf5_file, input_file, output_file)

    def run(self):
        # Compute dos
        from dcore.sumkdft_post import SumkDFTDCorePost
        sk = SumkDFTDCorePost(hdf_file=self.model_hdf5_file, use_dft_blocks=False, h_field=0.0)
        setup_sk(sk, 'w', self.params)
        results = {}
        results['dos'], results['dosproj'], results['dosproj_orb'] = \
            sk.dos_wannier_basis(broadening=self.params['broadening'],
                             mesh=self.params['mesh'],
                             with_Sigma=True, with_dc=self.params['with_dc'], save_to_file=False)

        if mpi.is_master_node():
            self.save_result(results)