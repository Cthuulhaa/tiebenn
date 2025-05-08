Velocity models
===============

This is a list of the velocity models currently implemented in TieBeNN. As further local, dedicated velocity models are published/discovered, they can be integrated in TieBeNN as specified further below.

List of available velocity models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - **Number assigned**
     - **Description**
   * - ``0``
     - IASP91 (`Kennett et al. 1995 <https://doi.org/10.1111/j.1365-246X.1995.tb03540.x>`_)
   * - ``1``
     - AK135 (`Montagner & Kennett 1996 <https://doi.org/10.1111/j.1365-246X.1996.tb06548.x>`_)
   * - ``2``
     - Layer-over-halfspace BGR velocity model (`Schlittenhardt 1999 <https://www.researchgate.net/profile/J-Schlittenhardt/publication/237600771_Regional_velocity_models_for_Germany_a_contribution_to_the_systematic_travel-time_calibration_of_the_international_monitoring_system/links/589dccbeaca272046aa92e2f/Regional-velocity-models-for-Germany-a-contribution-to-the-systematic-travel-time-calibration-of-the-international-monitoring-system.pdf>`_)
   * - ``3``
     - WET, velocity model for Germany and adjacent regions *(reference missing)*
   * - ``4``
     - WB2012, Vogtland/West Bohemia velocity model (1D, P+S), averaged from the 3D models of `Růžek & Horálek (2013) <https://doi.org/10.1093/gji/ggt295>`_
   * - ``5``
     - DEU, 3-layer model for Germany with a Moho depth of 28.5 km. A vp/vs ratio of 1.68 was used to generate S-wave velocities *(reference missing)*
   * - ``6``
     - Crust1.0 model based on epicenter location. Layers from `Leeman & Ammon <https://github.com/jrleeman/Crust1.0>`_
   * - ``7``
     - Crust1.0 + AK135: Crustal layers from Crust1.0, deeper structure from AK135
   * - ``8``
     - Insheim model (1D, P+S) from `Küperkoch et al. (2018) <https://doi.org/10.1785/0120170365>`_
   * - ``9``
     - Regional model for Landau-Insheim (by J. Borns)
   * - ``10``
     - Landau industry model (1D, P+S) for northern Upper Rhine Graben
   * - ``11``
     - 1D model for Central Alps (`Diehl et al. 2021 <https://doi.org/10.1029/2021JB022155>`_)
   * - ``12``
     - *(Deactivated)* 3D model for Central Alps (Diehl et al. 2021) – implementation under review
   * - ``13``
     - *(Deactivated)* 3D WEG model (shallow, high-res industry model around Rotenburg/Söhlingen/Soltau)
   * - ``14``
     - AlpsLocPS, 1D model for the Greater Alpine Region (`Brazsus et al. 2024 <https://doi.org/10.1093/gji/ggae077>`_)
   * - ``15``
     - BENS, 1D model for Northern Rhine Area (`Reamer & Hinzen 2004 <https://doi.org/10.1785/gssrl.75.6.713>`_)
   * - ``16``
     - ASZmod1, 1D model for Albstadt Shear Zone (`Mader et al. 2021 <https://doi.org/10.5194/se-12-1389-2021>`_)
   * - ``17``
     - *(Deactivated)* 3D P+S velocity model for the Upper Rhine Graben (`Lengliné et al. 2023 <https://doi.org/10.1093/gji/ggad255>`_)
   * - ``18``
     - PO1 (1D P+S) for Novy-Kostel swarm region, constrained (`Málek et al. 2023 <https://doi.org/10.1007/s00024-023-03250-w>`_)
   * - ``19``
     - PO1 + WB2012 (from 11 km depth)
   * - ``20``
     - KIT6 (1D P+S) for the Eastern Eifel Volcanic Field (`Ritter et al. 2024 <https://doi.org/10.1007/s10950-024-10257-w>`_)

When using ``--vel_mode manual``, the parameter ``--velmod`` must be set to one of the model numbers above.

When the **Crust1.0** model is selected, TieBeNN automatically extracts a 1D velocity profile at the epicenter coordinates.

Structure of NonLinLoc 1D velocity model files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each crustal velocity layer corresponds to one line in the text file using this format:

.. code-block:: text

   LAYER   DepthKm   VelP   GradVelP   VelS   GradVelS   Dens   GradDens

- ``VelP``, ``VelS`` in km/s
- ``GradVelP``, ``GradVelS`` in km/s/km (linear gradients)
- ``Dens`` in kg/m³, ``GradDens`` in kg/m³/km

.. note::

   The density is a carry-over field originally used in a waveform modeling code.
   It is not used in NonLinLoc programs, so any convenient value can be used
   *(A. Lomax, pers. comm.)*.

How to add a new 1D velocity model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Navigate to the folder: ``<repo>/tiebenn/data/velocity_models/``
   Create a new file named ``v_<number>`` (e.g. ``v21``), with one `LAYER` line per layer.

2. *(Optional)* List your model in the function docstring of ``velmods()`` in ``tools/velocity_models.py``.

3. Update/reinstall TieBeNN to integrate your new model (in the repository main directory):

   .. code-block:: bash

      pip install .

4. *(Optional / recommended)* Spread the love! Let the maintainer know via an issue or merge request so the model can be added to the official repo 🌍
