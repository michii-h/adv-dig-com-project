# Advanced Digital Communication Project (PStA)

The goal of the project is the implemetation of a _Software Defined Radio_ (SDR) based on the DRM standard, which works via a satellite communication link.
The SDR to be used is the ADLAM PLUTO (programed with MATLAB ) and the satellite to be used is the _Qatar-OSCAR 100_ (QO-100) (see: [amsat/qatar-oscar-100-qo-100](https://amsat-dl.org/eshail-2-amsat-phase-4-a-qatar-oscar-100/)).

The SDR has to fullfill the following **requirements**

- an extended simulation environment, which simulates the satellite communication channel (MATLAB Satellite Toolbox) with the appropriate QO-100 parameters, including the adjustment of the transmit power to Link Budget Estimation taking Antenna Gains and Noise margins into account
- Running the radio via the Adalm Pluto SDR
  - adopt the transmit parameters to the Narrowband Digital Mode of the QO-100 Transponder Bandplan (be careful with bandwidth and transmit frequencies; operation only with Mr. Lechner!!!) before being on air, check the transmit signal with the spectrum analyzer in the lab
  - Use a image of your choice as payload; include the call sign: DL0FHR
  - Optimize the fine synchronization with respect to SNR (FFT implementation)

The receiver shall be realized in Matlab (proposal: Live Script), which comprises a detailed (also mathematical) description of the different techniques and documentation of the corresponding code and figures (take care of labeling). The final presentation of the code shall be illustrated with a poster.

The grading depends on coding, functionality and the documentation!

- Submit the final Matlab Script "DRM_QO-100_YourName.m" or "DRM_QO-100_YourName.mlx" by uploading it to the Learning Campus
- Final presentation with Poster in the lab beginning of the summer semester 2025
- **Deadline** 15th of march 2025
