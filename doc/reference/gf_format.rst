File format for Green's function and self-energy
==================================================

The structure of the HDF5 dataset for Green's function and self-energy are as follows.


.. csv-table::
   :header: Data set, datatype, description
   :widths: 8, 10, 20

    __version, H5T_STRING, DCore_GfImFreq_v1
    data, "H5T_IEEE_F64LE (2*niw, N, N, 2)", "The first index runs from -niw to niw. The last one indices the real and imaginary parts"
    wn, "H5T_IEEE_F64LE (2*niw)",  "Imaginary frequenceis wn = (2*n+1)*pi*T for n = -niw, -niw+1, ..., niw-1"
