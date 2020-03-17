
***********
Run Modules
***********


pyPIC aims in implementing various simulation codes, whether they are Hybrid 
PIC of fully kinetic PIC codes. The following codes are implemented, they all
inherit from the base class Run. pyPIC implements the following codes: 

* :ref:`heckle`


:mod:`pypic.RunModules.run`
=============================

.. automodule:: pypic.RunModules.run
.. autoclass:: pypic.RunModules.run.Run


General information
-------------------

.. automethod:: pypic.RunModules.run.Run.display
.. automethod:: pypic.RunModules.run.Run.whereAmI



Getting physical quantities
---------------------------

Methods implemented in Run
^^^^^^^^^^^^^^^^^^^^^^^^^^

The following methods are general to all run object, they are therefore 
implemented in this base class and inherited by all run objects.


.. automethod:: pypic.RunModules.run.Run.GetVe
.. automethod:: pypic.RunModules.run.Run.GetVi
.. automethod:: pypic.RunModules.run.Run.GetVp
.. automethod:: pypic.RunModules.run.Run.GetNe
.. automethod:: pypic.RunModules.run.Run.GetNi
.. automethod:: pypic.RunModules.run.Run.GetNp
.. automethod:: pypic.RunModules.run.Run.GetP
.. automethod:: pypic.RunModules.run.Run.GetPe
.. automethod:: pypic.RunModules.run.Run.GetPp
.. automethod:: pypic.RunModules.run.Run.GetThermalEnergy



Methods to implement in run classes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following methods are not implemented in Run, however they must be 
implemented in all specific run object. pyPIC API will assume they are and 
other modules will possibly fail is they are not.


.. automethod:: pypic.RunModules.run.Run.GetV
.. automethod:: pypic.RunModules.run.Run.GetB
.. automethod:: pypic.RunModules.run.Run.GetE
.. automethod:: pypic.RunModules.run.Run.GetJ
.. automethod:: pypic.RunModules.run.Run.GetN
.. automethod:: pypic.RunModules.run.Run.GetFlux
.. automethod:: pypic.RunModules.run.Run.GetPxx
.. automethod:: pypic.RunModules.run.Run.GetPxy
.. automethod:: pypic.RunModules.run.Run.GetPxz
.. automethod:: pypic.RunModules.run.Run.GetPyy
.. automethod:: pypic.RunModules.run.Run.GetPyz
.. automethod:: pypic.RunModules.run.Run.GetPzz
.. automethod:: pypic.RunModules.run.Run.VgradV
.. automethod:: pypic.RunModules.run.Run.divP
.. automethod:: pypic.RunModules.run.Run.GetParticles
.. automethod:: pypic.RunModules.run.Run.GetMass
.. automethod:: pypic.RunModules.run.Run.GetCharge
.. automethod:: pypic.RunModules.run.Run.VgradV




.. _heckle:

:mod:`pypic.RunModules.heckle`
=============================

.. automodule:: pypic.RunModules.heckle
