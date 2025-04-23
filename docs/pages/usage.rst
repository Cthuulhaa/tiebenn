Usage
=====

This section shows an example which should make TieBeNN's usage clear.

Input file
~~~~~~~~~~

Input file syntax: ``YYYY-MM-DDTHH:MM:SS  latitude  longitude``

Example: ``2024-12-24T01:55:04 51.337  12.548``

Syntax
~~~~~~

With the created virtual environment activated, TieBeNN follows this syntax:

.. code-block:: bash

   tiebenn --event_file <EventFile> --max_epic_dist <MaxEpDist> --picker <Picker> --client <Client> --sds_dir <SDSDir> --min_detections <MinDetections> --plots <Plots> --vel_mode <VelMode> --velmod <VelMod> --ph_assoc <PhaseAssoc> --denoise <Denoise> --mult_windows <MultiWindows> --secs_before <SecsBef>

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - **Parameter**
     - **Description**
   * - ``--event_file``
     - Full path to input file with preliminary epicenter (latitude, longitude) and UTC datetime.
   * - ``--max_epic_dist``
     - Maximum epicentral distance (in km) for stations on which phase picks will be detected.
   * - ``--picker``
     - Select model for phase picking.  
       **EQTransformer**: ``sb_eqt``, ``seisbench_eqt``, ``seisbench_eqtransformer``, ``sb_eqtransformer``  
       **PhaseNet**: ``sb_pn``, ``seisbench_pn``, ``sb_phasenet``, ``seisbench_phasenet``  
       *(Not case sensitive)*.
   * - ``--client``
     - If set to ``SDS``, it accesses a SeisComp3 directory and fetches waveforms for stations within ``--max_epic_dist`` km.  
       Then attempts FDSN query. If set to ``FDSN``, it skips SDS lookup.
   * - ``--sds_dir``
     - Full path to SeisComp3 directory. Required if ``--client`` is set to ``SDS``.
   * - ``--min_detections``
     - Minimum number of stations with P- or S- picks required to complete the detection loop.
   * - ``--plots``
     - If ``True``, it plots waveforms with detections, phase picks sorted by distance, station maps,  
       and the confidence ellipsoid of the located event.
   * - ``--vel_mode``
     - Velocity model selection mode:  
       ``automatic`` (or ``auto``, ``a``): Picks model based on epicenter if available, or uses default  
       ``manual`` (or ``man``, ``m``): You must also define ``--velmod``.
   * - ``--velmod``
     - Number corresponding to the velocity model used in NonLinLoc.  
       :doc:`See the current velocity model list here <velocity_models>`.  
       *Note*: 3D models are supported, but not yet optimal for real-time use.
   * - ``--ph_assoc``
     - Phase associator: ``PyOcto`` (``pyocto``, ``p``) or ``GaMMA`` (``gamma``, ``g``). *(Not case sensitive)*
   * - ``--denoise``
     - If ``True``, applies DeepDenoiser to waveforms within 100 km of the epicenter.
   * - ``--mult_windows``
     - If ``True``, the phase picker searches in multiple time windows to improve ML stability.
   * - ``--secs_before``
     - If ``--mult_windows`` is ``False``, this sets how many seconds before the event the waveform starts.  
       Default is 0 seconds.

.. note::

   Commands related to the use of 3D NonLinLoc velocity models is deactivated, but they are still present in the code in case someone wants to pick its development up.

To locate the event with example UTC datetime and coordinates specified above, we type in the terminal (remember to activate the virtual environment created above):

.. code-block:: bash

   tiebenn --event_file <full_path_to_example_event> --max_epic_dist 150 --picker SeisBench_PhaseNet --client FDSN --min_detections 3 --plots True --vel_mode auto --ph_assoc PyOcto --denoise True --mult_windows True

Output
~~~~~~

If the location was successful, you’ll see something like:

.. code-block:: text

   ##############################################################
   NonLinLoc: Location completed.
   Origin time: 24-12-2024_01:55:4.517789
                Lat:51.328045 Long:12.543796 Depth: 22.685547 km
   --------------------------------------------------------------
   Location quality:
   RMS:0.285177 Number of phases: 130 Gap: 35.032 Distance from hypocenter to nearest station: 12.077 km
   Ellipsoid semi-major axis: 1.030196e+00
   ##############################################################

With ``--plots True``, the following figures are also generated:

* Phase picks on waveforms:

  .. image:: ../_static/example_picks.svg
     :alt: Phase picks
     :width: 500px
     :align: center

* Associated phases:

  .. image:: ../_static/example_phassoc.svg
     :alt: Phase association
     :width: 500px
     :align: center

* Map: epicenter and stations:

  .. image:: ../_static/example_epicenter_map.png
     :alt: Epicenter map
     :width: 500px
     :align: center

* Location confidence ellipsoid:

  .. image:: ../_static/example_ellipsoid.png
     :alt: Confidence ellipsoid
     :width: 500px
     :align: center

* Location Quality Score:

  .. image:: ../_static/example_LQS.svg
     :alt: LQS
     :width: 500px
     :align: center
