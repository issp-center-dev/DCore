import numpy
import time
import os
from triqs.gf import *
from dcore.tools import make_block_gf
from triqs_dft_tools import SumkDFT  # original SumkDFT
from dcore.sumkdft_opt import SumkDFT_opt  # optimized SumkDFT

def test_sumkdft_opt(request):
    org_dir = os.getcwd()
    os.chdir(request.fspath.dirname)

    kwargs = {
        'hdf_file': 'nis.h5',  # copied from test/pre_wannier
        'use_dft_blocks' : False,
        'h_field': 0.0,
    }
    
    print("Initialize original SumkDFT")
    sk_org = SumkDFT(**kwargs)
    print("Initialize optimized SumkDFT")
    sk_opt = SumkDFT_opt(**kwargs)
    print("Initialization done")
    
    n_corr_shells=sk_org.n_corr_shells
    print("n_corr_shells = {}".format(n_corr_shells))
    gf_struct = sk_org.gf_struct_solver
    print("gf_struct = {}".format(gf_struct))
    
    # set up Sigma_imp
    beta = 10
    n_iw = 100
    Sigma_imp = [make_block_gf(GfImFreq, gf_struct[icrsh], beta, n_iw) for icrsh in range(n_corr_shells)]
    numpy.random.seed(0)  # fix the seed of random number
    for icrsh in range(n_corr_shells):
        for sp, sigma in Sigma_imp[icrsh]:
            # Sigma is not hermitian but it doesn't matter
            sigma.data[:, :, :] = numpy.random.randn(*(sigma.data.shape)) \
                                + numpy.random.randn(*(sigma.data.shape)) * 1.0j
    
    # compute G_loc
    print("Compute G_loc")
    g_loc = []
    for sk, label in zip([sk_org, sk_opt], ['original ', 'optimized']):
        start = time.time()
        sk.set_Sigma(Sigma_imp)
        # sk.set_dc(dc_imp, dc_energ)
        # sk.set_mu(mu)
        g_loc.append(sk.extract_G_loc(with_Sigma=True, with_dc=False))
        print("{} : {} sec".format(label, time.time()-start))
    g_org = g_loc[0]
    g_opt = g_loc[1]
    
    # compare
    for icrsh in range(n_corr_shells):
        for sp in list(gf_struct[icrsh].keys()):
            assert numpy.allclose(g_org[icrsh][sp].data, g_opt[icrsh][sp].data)

    os.chdir(org_dir)