# NanOpticS
In-depth analysis of NANomaterials for OPTICal localized surface plasmon resonance Sensing
## Problems Solved
This software solves the problem of analysing the LSPR bands of plasmonic nanomaterials in monitoring processes, with hundreds of transmittance spectra, in a few seconds. It was conceived mainly for high-resolution LSPR spectroscopy systems. A statistical analysis of the LSPR bands with normalized spectral distributions is performed to compare different sensing platforms, and the signal-to-noise ratio (SNR) is calculated for each analysed parameter. Furthermore, all the spectra are fitted with a polynomial function, which enables a fast and direct analysis of the transmittance (LSPR band) minimum.
## Files
The instalation file - NANOPTICS24.exe

A text file containning the raw data to try the software - raw_transmittance_data.zip

This Readme file
## Installation
NANOPTICS must be installed under a 64-bit platform running Windows. 

During installation the user will be prompted to download (free) and install MATLAB Runtime R2018a (v 9.4).
## Running the software
To run the software, a file containing the transmittance spectra acquired over time, must have the wavelength (in nm) in the first column and the transmittance (in %), of the consecutive spectra, ordered in the next columns.

Some initial parameters must be input in the GUI:

1-The user starts by writing the "Monitoring name" that will be added to the produced figures and files to export;

2-The user must choose whether a "LSPR band" exists, or not;

3-If the LSPR band exists and is clearly seen, the user can choose to let the algorithm do an "Automatic wavelength range". If there is no clear observation of the LSPR band, or if the user prefers to choose another range, the wavelength range can be made manually;

4-The "Wavelength lower limit" should be set to include the relative maximum transmittance, while "Wavelength higher limit" should be set, when possible, as having the same y-axis coordinate of the wavelength lower limit;

5-The user may choose to "Save all spectra figures" for all the spectra analysed. It will show all transmittance spectra acquired during monitoring, as well as the analysed range with the correspondent fitting obtained using a polynomial function. If this option is not selected, only the first two spectra, along with last one, will be saved as figures. By default, this checkbox is not selected since this operation is time consuming;

6-The elapsed "Time between each spectrum" must be written, in milliseconds;

7-The user can choose to apply a linear "Drift correction" of the monitored parameters;

8- The total "Number of cycles" (each cycle with the reference and the test conditions) must be input. The last cycle of the sensing test must be done only with the reference condition, as it was programmed to be the “control” cycle. The time for each condition (half-cycle) should be the same, and the total cycle time should also be constant;

9-The user must enter the portion (percentage) of each "Half cycle analysis" of the monitoring to be considered for the average parameter shift calculation;

10-If the user is running a monitoring process to assess the refractive index sensitivity (RIS) of the sensor, the "Sensitivity" check box needs to be ticked, and then the refractive indices of the reference "P1" and test "P2" materials have to be input. 

During the analysis, a folder named “figures” and another one named “results” are created. The figures containing the individual spectra are saved in the first folder, and the figures containing all the calculated parameters in the second folder. All the files with the analysis results are also saved in the “results” folder.
## Version
2020a (2.4)
## License
This project is licensed under the MIT License
