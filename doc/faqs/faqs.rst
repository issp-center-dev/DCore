Frequently-Asked Questions
==========================

wien2k: FERMI ERROR when running `x lapw2 -almd -band`
------------------------------------------------------

In some versions of Wien2k, there is a problem in running `x lapw2 -almd -band`.

A hack solution is as follows:
1) `x lapw1 -band`
2) edit in2 file: replace 'TOT' with 'QTL', 'TETRA' with 'ROOT'
3) `x lapw2 -almd -band`
4) `dmftproj -band` (add the Fermi energy to file, it can be found by running `grep :FER *.scf`)

How do I plot the output of `spaghettis`?
-----------------------------------------

In python, you can do the following for example. You should pass the name of
the file written out by the spaghettis function.  Of course, you should change
the parameters as desired.

.. literalinclude:: plotting_spaghettis.py

x optic does not write a case.pmat file
---------------------------------------

Make sure that you set line 6 to "ON" and put a "1" to the following line.
The "1" is undocumented in Wien2k, but needed to have `case.pmat` written.
However, we are working on reading directly the `case.mommat2` file.

No module named pytriqs.*** error when running a script
-------------------------------------------------------

Make sure that have properly build, tested and installed TRIQS and DFTTools
using, make, make test and make install. Additionally, you should always
use pytriqs to call your scripts, e.g. pytriqs yourscript.py

Why is my calculation not working?
----------------------------------

Are you running in the right shell?
