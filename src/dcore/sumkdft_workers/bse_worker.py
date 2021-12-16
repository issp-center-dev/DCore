import numpy
import triqs.utility.mpi as mpi
from h5 import HDFArchive
from .worker_base import SumkDFTWorkerBase, setup_sk

class SumkDFTWorkerBSE(SumkDFTWorkerBase):
    """For computing BSE"""
    def __init__(self, model_hdf5_file, input_file, output_file) -> None:
        super().__init__(model_hdf5_file, input_file, output_file)

    def run(self):
        # chi0
        from dft_tools.sumk_dft_chi import SumkDFTChi
        # save div data (overwrite if data exist)
        if mpi.is_master_node():
            with HDFArchive(self.model_hdf5_file, 'a') as ar:
                if 'dft_input_chi' in ar:
                   del ar['dft_input_chi']
                ar.create_group('dft_input_chi')
                ar['dft_input_chi']['div'] = numpy.array(self.params['div'])
        # check if IBZ and FBZ data are saved separately
        dft_data_fbz = 'dft_input'
        if mpi.is_master_node():
            with HDFArchive(self.model_hdf5_file, 'r') as ar:
                if 'dft_input_fbz' in ar:
                    dft_data_fbz = 'dft_input_fbz'
        dft_data_fbz = mpi.bcast(dft_data_fbz)
        sk = SumkDFTChi(hdf_file=self.model_hdf5_file, use_dft_blocks=False, h_field=0.0,
                        dft_data_fbz=dft_data_fbz)
        setup_sk(sk, 'iwn', self.params)

        temp_file = None
        if self.params['use_temp_file']:
            temp_file = 'G_k_iw_temp.h5'

        sk.save_X0q_for_bse(list_wb=self.params['list_wb'],
                            n_wf_cutoff=self.params['n_wf_G2'],
                            qpoints_saved=self.params['X0q_qpoints_saved'],
                            h5_file=self.params['bse_h5_out_file'],
                            temp_file=temp_file,
                            nonlocal_order_parameter=False)