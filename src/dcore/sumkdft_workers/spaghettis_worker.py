from dcore._dispatcher import mpi
from .worker_base import SumkDFTWorkerBase, setup_sk

class SumkDFTWorkerSpaghettis(SumkDFTWorkerBase):
    """For computing Spaghettis"""
    def __init__(self, model_hdf5_file, input_file, output_file) -> None:
        super().__init__(model_hdf5_file, input_file, output_file)

    def run(self):
        from dcore.sumkdft_post import SumkDFTDCorePost
        sk = SumkDFTDCorePost(hdf_file=self.model_hdf5_file, use_dft_blocks=False, h_field=0.0,
                              bands_data=self.params['bands_data'])
        setup_sk(sk, 'w', self.params)
        results = {}
        results['akw'] = sk.spaghettis(broadening=self.params['broadening'],
            plot_range=None, ishell=None, save_to_file=None)

        if mpi.is_master_node():
            self.save_result(results)