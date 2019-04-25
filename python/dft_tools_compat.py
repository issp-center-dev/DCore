from pytriqs import version

major_version = int(version.version.split('.')[0])

if major_version == 1:
    from pytriqs.applications.dft.sumk_dft import SumkDFT, SumkDFTTools
elif major_version == 2:
    from triqs_dft_tools.sumk_dft import SumkDFT
    from triqs_dft_tools.sumk_dft_tools import SumkDFTTools
else:
    raise RuntimeError("Unknown triqs major version!")
