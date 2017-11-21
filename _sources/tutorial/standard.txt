Tutorial with preset model
==========================

t2g Bethe lattice
-----------------

:download:`dmft.ini <bethe-t2g/dmft.ini>`

.. literalinclude:: bethe-t2g/dmft.ini
                              
``pydmft_pre``
--------------

.. code-block :: bash

   $ pydmf_pre model.in

``pydmft``
----------

.. code-block :: bash

   $ pydmf dmft.ini bethe

``pydmft_post``
---------------

.. code-block :: bash

   $ pydmf_post dmft.ini bethe

