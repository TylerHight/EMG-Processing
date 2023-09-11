%% Main

% create an object to call the methods in signal_processing.m
sp = signal_processing;

fs = 1000; % sample frequency (Hz)

% file to be read
file_path = 'Output\Fatigue.xlsx';

% gets sample times and associated values
[st_cc, sig_cc] = sp.get_signal(file_path); 
st_cc = st_cc(1:length(sig_cc));

% generate filtered signal
sig_cc = sp.bandpass_butterworth(sig_cc, 4, [2 499.99], fs);

% 1) isolate the EMG data
%(sig, start time of interest, end time of interest, sample frequency)
[st, sig_cc] = sp.cut_sig_time(st_cc, sig_cc, 1, 4, fs);

name = "Fatigue";
figure("Name", name)
subplot(2, 1, 1)
plot(st, sig_cc)
title(name);
xlabel('time (s)')
ylabel('Volts')


% 2) Averaged Fourier power spectrum with using record lengths of 512 data points

% finds the number of 512 length sections that can be isolated
num_rows = floor(length(sig_cc)/512);
% shorten the array length into the largest possible multiple of 512
[st, sig] = sp.cut_sig_index(st, sig_cc, 1, 512*num_rows);
% reshape arrays into num_sections rows and 512 columns
st = reshape(st, [num_rows, 512]);
sig = reshape(sig, [num_rows, 512]);

% apply Fast Fourier Transform
[f_axis, fourier_t] = sp.fast_fourier_t(sig, fs);

% calculate Fourier power spectrum
fourier_t = mean(abs(fourier_t.^2));

% plot power spectrum
subplot(2, 1, 2);
plot(f_axis, fourier_t)
title("Fourier Power Spectrum with 512 points per window")
xlabel('Frequency (Hz)')
ylabel("Volts/Hz")


%% 3) Averaged spectra of 512 data point windows for each 10 second window of signal
% file to be read
file_path = 'Output\Fatigue.xlsx';

% gets sample times and associated values
[st, sig] = sp.get_signal(file_path); 

% generate filtered signal
sig = sp.bandpass_butterworth(sig, 4, [2 499.99], fs);

% number of data points in 10 seconds
np = 10*fs;

% finds the number of np length sections that can be isolated
num_rows = floor(length(sig)/np);

% shorten the array length into the largest possible multiple of np
[st, sig] = sp.cut_sig_index(st, sig, 1, np*num_rows);

% reshape arrays into num_sections rows and np columns
st = reshape(st, [np, num_rows]);
sig = reshape(sig, [np, num_rows]);
[f_axis, fourier_t] = sp.fast_fourier_t(sig, fs);
figure(5)
plot(f_axis, fourier_t)

% apply Fast Fourier Transform
[f_axis, fourier_t1] = sp.fast_fourier_t(sig, fs);

% calculate Fourier spectrum
fourier_s = mean(fourier_t1, 2);

% plot power spectrum
figure("Name", "Averaged Fourier Spectrum with 10 seconds per window")
sp.plot_fft(f_axis, fourier_s)


% 4) Mean frequency for each spectrum used in the average

% use only positive frequencies
f_axis_pos = f_axis((length(f_axis)/2)+1:length(f_axis));
f_axis_pos = reshape(f_axis_pos, [length(f_axis_pos), 1]);
fourier_ps_pos = fourier_t1(length(f_axis_pos):2*length(f_axis_pos)-1, :);

mean_f = sum(fourier_ps_pos.*f_axis_pos)./sum(fourier_ps_pos);
figure(10)
plot(0:10:length(mean_f)*10-1, mean_f)
title("Fatigue - Mean Frequency for Every 10 Second Window")
xlabel('Time (s)')
ylabel("Frequency (Hz)")



%% 5) Power in signals as a function of time

path = 'Output\Fatigue.xlsx';
fs = 1000;

% gets sample times and associated values
[stc, sigc] = sp.get_signal(path); 

% generate filtered signal
filtered = sp.bandpass_butterworth(sigc, 4, [2 499.99], fs);

% 1) isolate the EMG data
%(sample times, sig, start time of interest, end time of interest, 
% sample frequency)
%[stc, sigc] = sp.cut_sig_time(st1, filtered, 1, 13, fs);

% number of data points in 0.1 seconds
np = 0.1*fs;

% finds the number of np length sections that can be isolated
num_sections = floor(length(sigc)/np);

% shorten the array length into the largest possible multiple of np
[st, sigc] = sp.cut_sig_index(stc, sigc, 1, np*num_sections);

% reshape arrays into num_sections columns and np rows
sig = reshape(sigc, [np, num_sections]);

% find power in each 0.1s window
power = bandpower(sig);

% time array for 0.1 increments
t = 0:0.1:length(power)*0.1-(0.1);

figure(5)
plot(t, power)
title("Power in Fatigue - 0.1 second windows")
xlabel('Time (s)')
ylabel("Volts/Hz")


%% 6) determines whether the weight of the buckets affects the avg signal... 
% power and the mean frequency from the spectrum

sp = signal_processing;
fs = 1000;

% full bucket signal file path
file_path = 'Output\Force Effect (Isometric)\Full.xlsx';

% sample times and associated signal values
[isometric_ft1, isometric_f1] = sp.get_signal(file_path); 

% filter the signal
f_isometric_f1 = sp.bandpass_butterworth(isometric_f1, 4, [2 499.99], fs);

% isolate the EMG data by selecting start and end time
%(sig, start time of interest, end time of interest, sample frequency)
[st1, sig1] = sp.cut_sig_time(isometric_ft1, f_isometric_f1, 1, 14.5, fs);


% half full bucket signal file path
file_path = 'Output\Force Effect (Isometric)\Half.xlsx';

% sample times and associated signal values
[isometric_ft2, isometric_f2] = sp.get_signal(file_path); 

% filter the signal
f_isometric_f2 = sp.bandpass_butterworth(isometric_f2, 4, [2 499.99], fs);

% isolate the EMG data by selecting start and end time
%(sig, start time of interest, end time of interest, sample frequency)
[st, sig2] = sp.cut_sig_time(isometric_ft2, f_isometric_f2, 1, 14.5, fs);


% Determine whether the power signals have statistically significant 
% differences in average power. Returns 1 if different and 0 if not enough
% information to determine this.
diff = sp.compare_avg_power(sig1, sig2);


% Determining if the mean frequencies are significantly different:
dmf = sp.compare_mean_frequencies(sig1, sig2, fs);


%% 7) 
% determines if the isometric or concentric data differ in the avg signal 
% power and the mean frequency from the spectrum

sp = signal_processing;
fs = 1000;

% full bucket signal file path
path_full = 'Output\Force Effect (Isometric)\Full.xlsx';

% sample times and associated signal values
[isometric_ft1, isometric_f1] = sp.get_signal(file_path); 

% filter the signal
f_isometric_f1 = sp.bandpass_butterworth(isometric_f1, 4, [2 499.99], fs);

% isolate the EMG data by selecting start and end time
%(sig, start time of interest, end time of interest, sample frequency)
[st1, sig1] = sp.cut_sig_time(isometric_ft1, f_isometric_f1, 1, 14.5, fs);


% concentric signal file path
file_path = 'Output\Concentric Contraction.xlsx';
sp = signal_processing;
fs = 1000;

% sample times and associated signal values
[isometric_ft2, isometric_f2] = sp.get_signal(file_path); 

% filter the signal
f_isometric_f2 = sp.bandpass_butterworth(isometric_f2, 4, [2 499.99], fs);

% isolate the EMG data by selecting start and end time
%(sig, start time of interest, end time of interest, sample frequency)
[st2, sig2] = sp.cut_sig_time(isometric_ft2, f_isometric_f2, 2, 13, fs);


% Determine whether the power signals have statistically significant 
% differences in average power. Returns 1 if different and 0 if not enough
% information to determine this.
diff = sp.compare_avg_power(sig1, sig2);

% Determining if the mean frequencies are significantly different:
d = sp.compare_mean_frequencies(sig1, sig2, fs);


%% 8)
% Determine whether power changes with time  for the fatigue data by using 
% the p-value from a linear regression

file_path = 'Output\Fatigue.xlsx';
sp = signal_processing;
fs = 1000;

% sample times and associated signal values
[isometric_ft2, isometric_f2] = sp.get_signal(file_path); 

% filter the signal
f_isometric_f2 = sp.bandpass_butterworth(isometric_f2, 4, [2 499.99], fs);

sig = f_isometric_f2;
sig_power = sig.^2;
time = isometric_ft2;

% get p-value from linear regression
fit = fitlm(sig_power, time);
p_value = fit.Coefficients(2, 4);


% 9)
% Determine whether mean frequency changes with time for the fatigue data 
% by using the p-value from a linear regression using windows that increase
% in intervals of 10 seconds

sig_fatigue = f_isometric_f2;
time_fatigue = isometric_ft2;

% number of data points in 10 seconds
np = 10*fs;

% finds the number of np length sections that can be isolated
num_columns = floor(length(sig_fatigue)/np);

% shorten the array length into the largest possible multiple of np
[st1, sig1] = sp.cut_sig_index(time_fatigue, sig_fatigue, 1, np*num_columns);

% reshape arrays into np rows and num_columns columns (10 seconds of data
% per column)
st = reshape(st1, [np, num_columns]);
sig = reshape(sig1, [np, num_columns]);

% change the time intervals to every 10 seconds
st = st(1:np:end);

% calculate the mean frequency of each column
mean_f = meanfreq(sig, fs);

figure(9)
plot(st, mean_f)

% get p-value from linear regression
fit = fitlm(mean_f, st);
p_value = fit.Coefficients(2, 4);




