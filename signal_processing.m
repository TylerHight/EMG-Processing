% signal processing tools
classdef signal_processing
    methods
        function [t, sig] = EMG_signal(obj, file_path, fs)
            % sample times and associated signal values
            [t, sig] = obj.get_signal(file_path); 
            % filter the signal
            sig = obj.bandpass_butterworth(sig, 4, [2 499.99], fs);
        end


        function [reject] = compare_avg_power(~, sig1, sig2)
        % returns 0 if fail to reject and 1 if reject, meaning 1 if signals
        % have different means and 0 if not enough information to prove it

        % square the signals to get their power at any point
        sig1p = sig1.^2;
        sig2p = sig2.^2;

        reject = ttest2(sig1p, sig2p);
        end
        

        function [reject] = compare_mean_frequencies(obj, sig1, sig2, fs)
            % fft of signals
            [f1, fft1] = obj.fast_fourier_t(sig1, fs);
            [f2, fft2] = obj.fast_fourier_t(sig2, fs);
            
            % isolating positive frequencies %
            % sig1
            f1 = f1(f1>=0);
            f1 = reshape(f1, [length(f1), 1]);
            fft1 = fft1(length(f1):2*length(f1)-1, :);
            % sig2
            f2 = f2(f2>=0);
            f2 = reshape(f2, [length(f2), 1]);
            fft2 = fft2(length(f2):2*length(f2)-1, :);
            
            % Finding the mean of these arrays results in the mean frequencies. This 
            % expression differs from the expression used previously to find the mean
            % frequency in that the numerator is not summed and the denominator is
            % divided by the length of the arrays because these operations will be
            % performed in the ttest2() function. 
            in1 = (fft1.*f1)./(sum(fft1)./length(f1));
            in2 = (fft2.*f2)./(sum(fft2)./length(f2));
            
            % This determines if the mean frequencies are significantly different.
            % Since ttest2() calculates the mean of each array, and the mean of the 
            % arrays is the mean frequency, the 2-sample t-test is being performed on
            % the mean frequencies.
            reject = ttest2(in1, in2);
        end


        function [sample_times, signal] = get_signal(~, file_path)
            % get the signal values
            signal = xlsread(file_path, 'B:B');
            % get the time of each signal recording
            sample_times = xlsread(file_path, 'A:A');
        end
        
        
        function [filtered] = bandpass_butterworth(~, signal, order, f_range, fs)
            % (signal outputs, order of filter, [low_f high_f], sample frequency)
            % first parameter: order 2n filter
            % second paramter: cutoff frequencies
            % butter() creates 2nth order filter where n is the parameter
            input_o = order/2; 
            % transfer function coefficients
            [b, a] = butter(input_o, f_range/(fs/2), "bandpass"); 
            filtered = filter(b, a, signal); % apply the Butterworth filter
        end
        
        
        function [st, sig] = cut_sig_time(~, sample_times, signal, start_time, end_time, sample_f)
            % cuts out a part of the signal based on time
        
            % index of the start of the wanted part of signal
            i_start = sample_f*start_time;
            % index of the end of the wanted part of signal
            i_end = sample_f*end_time;
            st = sample_times(i_start:i_end);
            sig = signal(i_start:i_end);
        end
        
        
        function [st, sig] = cut_sig_index(~, sample_times, signal, i_start, i_end)
            % cuts out a part of the signal based on index
            st = sample_times(i_start:i_end);
            sig = signal(i_start:i_end);
        end
        
        
        function [fs] = find_sample_frequency(~, sample_times)
            % finds the sample frequency for an input sample time array
        
            tbs = mean(diff(sample_times)); % average time between measurements
            fs = round(1/tbs); % sample frequency
        end
        
        
        function [faxis, ft_shift] = fast_fourier_t(~, signal_array, sample_freq)
            %%% generate FFT for sampled signal
            num_samples = length(signal_array);
            % generate fast Fourier transform (FFT) and circshift
            fast_f_transform = fft(signal_array)/num_samples;
            % magnitude axis
            ft_shift = circshift(fast_f_transform, round(num_samples/2));
            % generate frequency axis (faxis)
            faxis = -sample_freq/2:sample_freq/num_samples:sample_freq/2-...
                sample_freq/num_samples;
        end
        
        
        function [] = plot_fft(~, faxis, ftshift)
            % plots the abs of the fft with a linear x-axis and on a log-log plot
            subplot(2, 1, 1);
            plot(faxis, abs(ftshift))
            xlabel('Frequency (Hz)');
            ylabel('Magnitude');
            title("Spectra Averaged Across 10 Second Windows")
            subplot(2, 1, 2);
            loglog(faxis, abs(ftshift))
            title("Spectra Averaged Across 10 Second Windows - log-log")
            xlabel('Frequency (Hz)');
            ylabel('Magnitude');
        end

        
        function [] = plot_power_spectrum(~, faxis, ftshift)
            % plots the abs of the fft with a linear x-axis and on a log-log plot
            subplot(2, 1, 1);
            plot(faxis, abs(ftshift))
            xlabel('Frequency (Hz)');
            ylabel('(V^2)/Hz');
            subplot(2, 1, 2);
            loglog(faxis, abs(ftshift))
            xlabel('Frequency (Hz)');
            ylabel('(V^2)/Hz');
        end

    end
end