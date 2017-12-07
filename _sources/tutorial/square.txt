Tutorial with single-band 2D Hubbard model
==========================================

:download:`dmft_square.ini <square/dmft_square.ini>`

.. literalinclude:: square/dmft_square.ini
                              
``dcore_pre``
--------------

.. code-block :: bash

   $ pydmf_pre dmft_square.ini

``dcore``
----------

.. code-block :: bash

   $ pydmf dmft_square.ini

``dcore_post``
---------------

.. code-block :: bash

   $ pydmf_post dmft_square.ini
   $ gnuplot square_akw.gp

.. image:: square/akw.png
   :width: 700
   :align: center

.. code-block :: gnuplot

   gnuplot> set xlabel "Energy"
   gnuplot> set ylabel "DOS"
   gnuplot> plot "square_dos.dat" w l

.. image:: square/dos.png
   :width: 700
   :align: center
