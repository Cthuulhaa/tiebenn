Workflow
========

Using the coordinates (latitude and longitude) and the UTC time of a local event as input, TieBeNN processes seismic event locations through the following stages:

1. **Waveform data fetching**: A catalog of stations around the epicenter is produced. Then, waveform data in miniSEED format are retrieved using an ObsPy client on FDSN servers or from a SDS directory structure.

2. **Waveform data preprocessing**: Empty channels are removed, masked channels are split, and data are detrended and bandpass filtered. Optionally, stations within 100 km from the epicenter are denoised using the DeepDenoiser model `(Zhu et al. 2019) <https://arxiv.org/abs/1811.02695>`_.

3. **Phase picking**: P- and S-phases (first arrivals) are detected using phase-picking models (either EQTransformer or PhaseNet, available in the `SeisBench toolbox <https://github.com/seisbench/seisbench>`_.

4. **Phase association**: Detected picks are passed through phase associators to discard possible false detections. TieBeNN supports `PyOcto <https://github.com/yetinam/pyocto>`_ and `GaMMA <https://github.com/AI4EPS/GaMMA>`_.

5. **Export outputs**: Event-associated phase picks are exported to CSV files, including station coordinates, arrival times, pick probabilities, and signal-to-noise ratios. NonLinLoc input files are generated. Optionally, waveform and pick figures per station, as well as phase association plots, are produced.

6. **Probabilistic hypocenter estimation**: The generated files are used by NonLinLoc for hypocenter estimation. Optionally, the following figures are created:

  i. Epicenter and stations on a map
  ii. Waveforms with picks sorted by epicentral distance
  iii. location confidence ellipsoid.

7. **Location quality assessment**: Location metrics are gathered to compute the Location Quality Score (LQS) and to generate a visualization. **A description of this metric should be available in a manuscript, currently in preparation**

.. note::

  TieBeNN loops through this process until the minimum required detections within a given epicentral distance are obtained. If not, the search radius is gradually expanded. If phase picks are insufficient within 200 km, the run ends with an unsuccessful event location.
