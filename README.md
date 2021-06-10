# FIR_Hamming_Filter
FIR filter based on a Hamming window

C# Library to calculate the coefficients of a FIR filter with Hamming window and to filter the data

This code calculates the coefficients of the Band-pass, Band-stop, Low-pass and High-pass FIR filters. It also filters the data, but no zero-phase delay is applied.

Each filter function will return a 2 rows x N coefficients 2D vector, where Row 1 = Numerator and Row 2 = Denumerator. 
Band-pass: the function is "double[][] Fir_LP (int filt_order, double sf, double f1, double f2)". The first argument if the order of the filter, the second one if the sampling frequency, while the last two arguments are the normalized two cut-off frequencies (f1/(sf/2), f2/(sf/2)), which means that the cutoff frequencies must be within the interval of (0,1);

Band-stop: the function is "double[][] Fir_LP (int filt_order, double sf, double f1, double f2)". The first argument if the order of the filter, the second one if the sampling frequency, while the last two arguments are the normalized two cut-off frequencies (f1/(sf/2), f2/(sf/2)), which means that the cutoff frequencies must be within the interval of (0,1). If the order of the filter is odd, it will be increased by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency;
Low-pass: the function is "double[][] Fir_LP(int filt_order, double sf, double f1)". The first argument if the order of the filter, the second one if the sampling frequency, while the last one is the normalized cut-off frequency (f1/(sf/2)), which means that the cutoff frequency must be within the interval of (0,1);
Low-pass: the function is "double[][] Fir_HP(int filt_order, double sf, double f1)". The first argument if the order of the filter, the second one if the sampling frequency, while the last one is the normalized cut-off frequency (f1/(sf/2)), which means that the cutoff frequency must be within the interval of (0,1). If the order of the filter is odd, it will be increased by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency;
Filter the data: the method is "double[] Filter_Data(double[][] coeff_filt, double[] pre_filt_signal)". The two arguments are the filter coefficients and the signal to be filtered. It returns the filtered signal.

If you have any question and/or want to report bugs, please e-mail me (Ale) at: pressalex@hotmail.com
