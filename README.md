# test2
# Author: Sean Maher
# Date: April 2020

If you find this repository useful in your own work please cite Maher et 
al., (2020).

Maher, S., Matoza, R., de Groot-Hedlin, C., Gee. K., Fee, D., Yokoo, A., 2020, 
Investigating Spectral Distortion of Local Volcano Infrasound by Nonlinear Propagation at 
Sakurajima Volcano, Japan: Journal of Geophysical Research, doi: 10.1029/2019JB018284.

--

This folder contains an example of the workflow involved in estimating nonlinear energy
transfer in high-amplitude acoustic waves using the frequency-domain Generalized Burgers 
Equation described by Reichman et al., (2016) as implemented by Maher et al. (2020) for 
volcano infrasound. The data provided come from a five-station deployment of infrasound 
sensors at Sakurajima Volcano in July 2013, as described by Fee et al., (2014). This 
example closely reproduces the results shown in Figure 3a-3d of Maher et al., (2020) 
including waveforms, a spectrogram, power spectral density curves, and the nu_N_tot curves. 

The folder "data" contains day-long infrasound records in MSEED format for five stations
at Sakurajima on 7/22/2013. The data can be freely obtained through IRIS DMC (network code 
YO).

The .npy files contain power spectra for the 5th, 50th, 90th, and 95th percentiles of noise 
averaged across the YO network as described by Maher et al., (2020).

The Python script 1_make_ascii.py generates text files containing waveforms for the example
event at each station (*.txt) and smoothed spectra for power spectral density and nu_N 
(*.dat). This script automatically triggers the plotting script 2_plot_results.py. 

Note that the resulting spectra do not precisely replicate the results in Figure 3a-3d of 
Maher et al., (2020). The spectra in the paper were created using a multi taper code written 
in Fortran. The spectra in this example are made using the Python library NiTime 
(http://nipy.org/nitime/index.html) and smoothed using a LOWESS function (Locally Weighted 
Scatterplot Smoothing) from the module statsmodels 
(https://www.statsmodels.org/stable/index.html). The similarity of the PSD and nu_N
spectra between the two methods suggests that the NiTime function is robust. The NiTime method
is slower than the Fortran code, but is used here for convenience in demonstrating the method 
using a fully Python workflow.



References:

Fee, D., Yokoo, A., & Johnson, J. B. (2014). Introduction to an open community infrasound dataset 
from the actively erupting Sakurajima Volcano, Japan. Seismological Research Letters, 85(6), 
1151–1162. https://doi.org/10.1785/0220140051

Maher, S., Matoza, R., de Groot-Hedlin, C., Gee. K., Fee, D., Yokoo, A., 2020 (in press), 
Investigating Spectral Distortion of Local Volcano Infrasound by Nonlinear Propagation at 
Sakurajima Volcano, Japan: Journal of Geophysical Research, doi: 10.1029/2019JB018284.

Reichman, B. O., Gee, K. L., Neilsen, T. B., & Miller, K. G. (2016). Quantitative analysis of a 
frequency-domain nonlinearity indicator. The Journal of the Acoustical Society of America, 139(5), 
2505–2513. https://doi.org/10.1121/1.4945787
