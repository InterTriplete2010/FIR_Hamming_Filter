using System;
using System.Numerics;
using System.Collections.Generic;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using System.Net.Http.Headers;
using System.Linq;
using MathNet.Numerics;

namespace FIR_HAMM
{
    public class FIR_HAMM
    {

        //Global variables
        int nbands;
        int first_band;
        double[] freq;
        double[] mags;
        double[] magnitude;
        double[] hamm_window;
        double[] F_Matlab;
        double[] h_Matlab;
        double[][] b_filter;


        //-------------------------------------------------------------------------------------------------------------------------------------//
        //Step 1: Compute the vector of frequencies to pass to FIRLS
        private void Desired_Freq(double Wnf_1)
        {

            nbands = 2;
            freq = new double[nbands * 2];

            //create the pair of frequency points based on Wn
            freq[0] = 0;
            freq[1] = Wnf_1;
            freq[2] = Wnf_1;
            freq[3] = 1;
        }

        //Overload of Step 1
        private void Desired_Freq(double Wnf_1, double Wnf_2)
        {

            nbands = 3;
            freq = new double[nbands * 2];

            //create the pair of frequency points based on Wn
            freq[0] = 0;
            freq[1] = Wnf_1;
            freq[2] = Wnf_1;
            freq[3] = Wnf_2;
            freq[4] = Wnf_2;
            freq[5] = 1;
        }
        //-------------------------------------------------------------------------------------------------------------------------------------//

        //-------------------------------------------------------------------------------------------------------------------------------------//
        //Step 2: Compute the magnitude vector
        private void Desired_MAG(int type_filt, int nbands)
        {

            if (type_filt == 1 || type_filt == 2)
            {

                first_band = 0;

            }

            else
            {

                first_band = 1;

            }


            mags = new double[nbands];

            for (int kk = 0; kk < nbands; kk++)
            {

                mags[kk] = (first_band + kk) % 2;

            }

            int track_mags = 0;
            magnitude = new double[nbands * 2];

            for (int kk = 0; kk < nbands * 2; kk += 2)
            {

                magnitude[kk] = mags[track_mags];
                magnitude[kk + 1] = mags[track_mags];

                track_mags++;

            }

        }

        //Step 3: Build a Hamming Window
        private void Hamm_win(int filter_order)
        {

            double coeff_I = 0.54;
            double coeff_II = 0.46;

            hamm_window = new double[filter_order];

            for (int kk = 0; kk < filter_order; kk++)
            {

                hamm_window[kk] = coeff_I - coeff_II * Math.Cos(2 * Math.PI * kk / (filter_order - 1));

            }

        }
        //-------------------------------------------------------------------------------------------------------------------------------------//

        //-------------------------------------------------------------------------------------------------------------------------------------//
        //Step 4: Sinc function
        private double Sinc_Function(double input_value)
        {
            double output_sinc;

            if (input_value == 0)
            {

                output_sinc = 1;

            }

            else
            {

                output_sinc = Math.Sin(input_value * Math.PI) / (input_value * Math.PI);

            }

            return output_sinc;

        }
        //-------------------------------------------------------------------------------------------------------------------------------------//

        //-------------------------------------------------------------------------------------------------------------------------------------//
        //Step 5: Compute windowed impulse response
        private void firls(int filter_order, int filter_type, double[] freq, double[] magnitude)
        {
            bool odd_order = false;

            if ((filter_order % 2) != 0)
            {

                filter_order += 1;
                odd_order = true;

            }

            double freq_length = freq.Length;
            double length_weight = Math.Floor(freq_length / 2);
            double[] weight = new double[Convert.ToInt32(length_weight)];
            for (int kk = 0; kk < weight.Length; kk++)
            {

                weight[kk] = 1;

            }

            F_Matlab = new double[freq.Length];
            for (int kk = 0; kk < F_Matlab.Length; kk++)
            {

                F_Matlab[kk] = freq[kk] / 2;

            }

            double[] amp_Matlab = new double[magnitude.Length];
            amp_Matlab = magnitude;

            double[] wt_Matlab = new double[Convert.ToInt32(length_weight)];
            for (int kk = 0; kk < wt_Matlab.Length; kk++)
            {

                wt_Matlab[kk] = Math.Abs(Math.Sqrt(weight[kk]));

            }

            //Calculate the "diff" of F_Matlab
            double[] dF_Matlab = new double[freq.Length - 1];
            for (int kk = 0; kk < dF_Matlab.Length; kk++)
            {

                dF_Matlab[kk] = F_Matlab[kk + 1] - F_Matlab[kk];

            }

            int lendF = freq.Length - 1;

            // validate weight
            double[] temp_W = new double[wt_Matlab.Length];
            for (int kk = 0; kk < temp_W.Length; kk++)
            {

                temp_W[kk] = wt_Matlab[kk] - wt_Matlab[0];

            }

            // find the order
            int L_Matlab = filter_order / 2;

            // odd order
            bool Nodd = false;

            if (odd_order == false)
            {

                Nodd = (((filter_order + 1) % 2) == 1);

            }

            else
            {

                Nodd = (((filter_order) % 2) == 1);

            }

            // initialize b0
            double b0_Matlab = 0;

            filter_type = 0;

            // Type I and Type II linear phase FIR
            double[] m_Matlab;
            if (odd_order == true)
            {

                m_Matlab = new double[L_Matlab];

            }

            else
            {

                m_Matlab = new double[L_Matlab + 1];

            }

            if (filter_type == 0)
            {
                // Basis vectors are cos(2 * pi * m * f)
                if (!Nodd == true)
                {

                    for (int kk = 0; kk < m_Matlab.Length; kk++)
                    {

                        m_Matlab[kk] = kk + 0.5;

                    }

                }

                else
                {

                    for (int kk = 0; kk < m_Matlab.Length; kk++)
                    {

                        m_Matlab[kk] = kk;

                    }

                }

            }

            double[] k_Matlab;
            if (odd_order == true)
            {

                k_Matlab = new double[m_Matlab.Length];

            }

            else
            {

                k_Matlab = new double[m_Matlab.Length - 1];

            }


            if (Nodd == true)
            {

                for (int kk = 0; kk < k_Matlab.Length; kk++)
                {

                    k_Matlab[kk] = m_Matlab[kk + 1];

                }

                b0_Matlab = 0; //First entry must be handled separately(where k(1) = 0)

            }

            else
            {

                k_Matlab = m_Matlab;

            }

            // preallocate b matrix
            double[] b_Matlab = new double[k_Matlab.Length];
            double m_s_Matlab;
            double b1_Matlab;
            double[] A_Matlab = new double[magnitude.Length];
            A_Matlab = amp_Matlab;

            for (int kk = 0; kk < F_Matlab.Length; kk += 2)
            {

                m_s_Matlab = (A_Matlab[kk + 1] - A_Matlab[kk]) / (F_Matlab[kk + 1] - F_Matlab[kk]);   // Slope
                b1_Matlab = A_Matlab[kk + 1] - m_s_Matlab * F_Matlab[kk];                           // y - intercept

                if (Nodd == true)
                {

                    b0_Matlab = b0_Matlab + (b1_Matlab * (F_Matlab[kk + 1] - F_Matlab[kk]) + m_s_Matlab / 2 * (F_Matlab[kk + 1] * F_Matlab[kk + 1] - F_Matlab[kk] * F_Matlab[kk])) * (Math.Pow((wt_Matlab[(kk + 2) / 2 - 1]), 2));

                }

                for (int ll = 0; ll < b_Matlab.Length; ll++)
                {

                    b_Matlab[ll] += (m_s_Matlab / (4 * Math.PI * Math.PI) * (Math.Cos(2 * Math.PI * k_Matlab[ll] * F_Matlab[kk + 1]) - Math.Cos(2 * Math.PI * k_Matlab[ll] * F_Matlab[kk])) / (k_Matlab[ll] * k_Matlab[ll])) * (Math.Pow(wt_Matlab[(kk + 2) / 2 - 1], 2));
                    b_Matlab[ll] += (F_Matlab[kk + 1] * (m_s_Matlab * F_Matlab[kk + 1] + b1_Matlab) * Sinc_Function(2 * k_Matlab[ll] * F_Matlab[kk + 1]) - F_Matlab[kk] * (m_s_Matlab * F_Matlab[kk] + b1_Matlab) * Sinc_Function(2 * k_Matlab[ll] * F_Matlab[kk])) * (Math.Pow(wt_Matlab[(kk + 2) / 2 - 1], 2));
                }

            }

            //Increase the size of the array and add b0_Matlab in the first cell if Nodd is true
            if (Nodd == true)
            {
                Array.Resize(ref b_Matlab, b_Matlab.Length + 1);

                for (int hh = 0; hh < b_Matlab.Length; hh++)
                {

                    if (hh == b_Matlab.Length - 1)
                    {

                        b_Matlab[b_Matlab.Length - hh - 1] = b0_Matlab;

                    }

                    else
                    {

                        b_Matlab[b_Matlab.Length - hh - 1] = b_Matlab[b_Matlab.Length - hh - 2];

                    }

                }
            }

            double[] a_Matlab = new double[b_Matlab.Length];
            for (int hh = 0; hh < a_Matlab.Length; hh++)
            {

                a_Matlab[hh] = Math.Pow(wt_Matlab[0], 2) * 4 * b_Matlab[hh];

                if (hh == 0 && Nodd == true)
                {

                    a_Matlab[hh] /= 2;

                }

            }

            //-------------------------------------------------------------------------------------//
            //Compute the unwindowed impulse response of the filter
            if (Nodd == true)
            {

                h_Matlab = new double[filter_order + 1];

                for (int hh = 0; hh < h_Matlab.Length; hh++)
                {

                    if (hh < filter_order / 2)
                    {

                        h_Matlab[hh] = a_Matlab[a_Matlab.Length - hh - 1] / 2;


                    }

                    if (hh == filter_order / 2)
                    {

                        h_Matlab[hh] = a_Matlab[0];

                    }

                    if (hh > filter_order / 2)
                    {

                        h_Matlab[hh] = a_Matlab[hh - h_Matlab.Length / 2] / 2;

                    }

                }

                //-------------------------------------------------------------------------------------//
                //Compute the windowed (Hamming) impulse response of the filter
                Hamm_win(h_Matlab.Length);
                for (int kk = 0; kk < h_Matlab.Length; kk++)
                {

                    h_Matlab[kk] = h_Matlab[kk] * hamm_window[kk];

                }
                //-------------------------------------------------------------------------------------//

            }

            else
            {

                if (odd_order == true)
                {

                    h_Matlab = new double[filter_order];

                }

                else
                {

                    h_Matlab = new double[filter_order + 1];

                }

                for (int hh = 0; hh < h_Matlab.Length; hh++)
                {

                    if (hh < h_Matlab.Length / 2)
                    {

                        h_Matlab[hh] = a_Matlab[a_Matlab.Length - hh - 1] * 0.5;

                    }

                    else
                    {

                        h_Matlab[hh] = a_Matlab[hh - a_Matlab.Length] * 0.5;

                    }

                }

                //-------------------------------------------------------------------------------------//
                //Compute the windowed (Hamming) impulse response of the filter
                Hamm_win(h_Matlab.Length);
                for (int kk = 0; kk < h_Matlab.Length; kk++)
                {

                    h_Matlab[kk] = h_Matlab[kk] * hamm_window[kk];

                }
                //-------------------------------------------------------------------------------------//

            }

        }

        //-------------------------------------------------------------------------------------//
        //Step 6: Scale the filter coefficients
        private void Scale_Filter()
        {
            double f0_Matlab;
            //b_filter = new double[h_Matlab.Length];
            b_filter = new double[2][];

            for (int ff = 0; ff < 2; ff++)
            {

                b_filter[ff] = new double[h_Matlab.Length];

            }

            //Setting the first value of the denumerator to "1"
            b_filter[1][0] = 1;

            if (first_band == 1)
            {

                double temp_sum = 0;

                for (int hh = 0; hh < h_Matlab.Length; hh++)
                {

                    temp_sum += h_Matlab[hh];

                }

                for (int hh = 0; hh < h_Matlab.Length; hh++)
                {

                    b_filter[0][hh] = h_Matlab[hh] / temp_sum;

                }

            }

            else
            {

                if (freq[3] == 1)
                {

                    // Unity gain at Fs / 2
                    f0_Matlab = 1;

                }

                else
                {

                    // unity gain at center of first passband
                    f0_Matlab = (freq[2] + freq[3]) / 2;

                }

                Complex complex_imag = new Complex(0.0, 1.0);
                Complex temp_den = new Complex(0.0, 0.0);
                for (int kk = 0; kk < h_Matlab.Length; kk++)
                {

                    temp_den += Complex.Exp(-complex_imag * 2 * Math.PI * kk * (f0_Matlab / 2)) * h_Matlab[kk];

                }

                double temp_den_double = Complex.Abs(temp_den);

                for (int kk = 0; kk < h_Matlab.Length; kk++)
                {

                    b_filter[0][kk] = h_Matlab[kk] / temp_den_double;

                }
            }
        }
        //-------------------------------------------------------------------------------------------------------------------------------------//

        //-------------------------------------------------------------------------------------------------------------------------------------//
        //Calculate the coefficients
        //Low-Pass filter
        public double[][] Fir_LP(int filt_order, double sf, double f1)
        {

            int filter_type = 0;
            //Normalizing the corner frequency with respect to the nyquist frequency
            f1 = f1 / (sf / 2);

            //Check that the normalized frequencies are within the correct range of values
            if (f1 <= 0 | f1 >= 1)
            {

                throw new Exception("Cut-off frequencies must be in the (0,1) range");

            }

            //Check that the order of the filter is > 0
            if (filt_order <= 0)
            {

                throw new Exception("The order of the filter must be > 0");

            }

            Desired_Freq(f1);
            Desired_MAG(filter_type, nbands);
            firls(filt_order, filter_type, freq, magnitude);
            Scale_Filter();

            return b_filter;

        }

        //High-Pass filter
        public double[][] Fir_HP(int filt_order, double sf, double f1)
        {

            int filter_type = 1;

            int temp_odd = filt_order % 2;
            if (temp_odd == 1)
            {

                //Increasing the order of the filter by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency
                filt_order += 1;

            }

            //Normalizing the corner frequency with respect to the nyquist frequency
            f1 = f1 / (sf / 2);

            //Check that the normalized frequencies are within the correct range of values
            if (f1 <= 0 | f1 >= 1)
            {

                throw new Exception("Cut-off frequencies must be in the (0,1) range");

            }

            //Check that the order of the filter is > 0
            if (filt_order <= 0)
            {

                throw new Exception("The order of the filter must be > 0");

            }

            Desired_Freq(f1);
            Desired_MAG(filter_type, nbands);
            firls(filt_order, filter_type, freq, magnitude);
            Scale_Filter();

            return b_filter;

        }

        //Band-Pass filter
        public double[][] Fir_BP(int filt_order, double sf, double f1, double f2)
        {

            int filter_type = 2;

            //Normalizing the corner frequencies with respect to the nyquist frequency
            f1 = f1 / (sf / 2);
            f2 = f2 / (sf / 2);

            //Check that the normalized frequencies are within the correct range of values
            if (f1 <= 0 | f1 >= 1 | f2 <= 0 | f2 >= 1)
            {

                throw new Exception("Cut-off frequencies must be in the (0,1) range");

            }

            //Check that f1 < f2
            if (f1 > f2)
            {

                throw new Exception("The corner frequency f1 needs to be lower than f2");

            }

            //Check that the order of the filter is > 0
            if (filt_order <= 0)
            {

                throw new Exception("The order of the filter must be > 0");

            }

            Desired_Freq(f1, f2);
            Desired_MAG(filter_type, nbands);
            firls(filt_order, filter_type, freq, magnitude);
            Scale_Filter();

            return b_filter;

        }

        //Band-Stop filter
        public double[][] Fir_BS(int filt_order, double sf, double f1, double f2)
        {

            int filter_type = 3;

            int temp_odd = filt_order % 2;
            if (temp_odd == 1)
            {

                //Increasing the order of the filter by 1, as odd order symmetric FIR filters must have a gain of zero at the Nyquist frequency
                filt_order += 1;

            }

            //Normalizing the corner frequencies with respect to the nyquist frequency
            f1 = f1 / (sf / 2);
            f2 = f2 / (sf / 2);

            //Check that the normalized frequencies are within the correct range of values
            if (f1 <= 0 | f1 >= 1 | f2 <= 0 | f2 >= 1)
            {

                throw new Exception("Cut-off frequencies must be in the (0,1) range");

            }

            //Check that f1 < f2
            if (f1 > f2)
            {

                throw new Exception("The corner frequency f1 needs to be lower than f2");

            }

            //Check that the order of the filter is > 0
            if (filt_order <= 0)
            {

                throw new Exception("The order of the filter must be > 0");

            }

            Desired_Freq(f1, f2);
            Desired_MAG(filter_type, nbands);
            firls(filt_order, filter_type, freq, magnitude);
            Scale_Filter();

            return b_filter;

        }
        //-------------------------------------------------------------------------------------------------------------------------------------//

        //-------------------------------------------------------------------------------------------------------------------------------------//
        //Filter the data by using the Direct-Form II Transpose, as explained in the Matlab documentation
        public double[] Filter_Data(double[][] coeff_filt, double[] pre_filt_signal)
        {

            double[] filt_signal = new double[pre_filt_signal.Length];
            Array.Clear(filt_signal, 0, filt_signal.Length);

            double[][] w_val = new double[coeff_filt[0].Length][];


            for (int ff = 0; ff < coeff_filt[0].Length; ff++)
            {

                w_val[ff] = new double[pre_filt_signal.Length];

            }


            //Convolution product to filter the data
            for (int kk = 0; kk < pre_filt_signal.Length; kk++)
            {

                if (kk == 0)
                {

                    filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0];


                    for (int ww = 1; ww < coeff_filt[0].Length; ww++)
                    {

                        w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];


                    }

                }

                else
                {

                    filt_signal[kk] = pre_filt_signal[kk] * coeff_filt[0][0] + w_val[0][kk - 1];

                    for (int ww = 1; ww < coeff_filt[0].Length; ww++)
                    {

                        w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] + w_val[ww][kk - 1] - filt_signal[kk] * coeff_filt[1][ww];

                        if (ww == coeff_filt[0].Length - 1)
                        {

                            w_val[ww - 1][kk] = pre_filt_signal[kk] * coeff_filt[0][ww] - filt_signal[kk] * coeff_filt[1][ww];

                        }

                    }

                }

            }

            return filt_signal;

        }
        //-------------------------------------------------------------------------------------------------------------------------------------//

    }
}

