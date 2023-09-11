% signal processing object
sp = signal_processing;
% sample freq (Hz)
fs = 1000;

% File paths
fp_iso_full = 'Output\Force Effect (Isometric)\Full.xlsx';
fp_iso_half = 'Output\Force Effect (Isometric)\Half.xlsx';
% maximum voluntary contraction data
mvc_1 = 'Output\Maximum Voluntary Contractions\1MVC';
mvc_2 = 'Output\Maximum Voluntary Contractions\2MVC';
mvc_3 = 'Output\Maximum Voluntary Contractions\3MVC';
cc = 'Output\Concentric Contraction.xlsx';
fat = 'Output\Fatigue.xlsx';

files = NaN(400);
files = [fp_iso_full; fp_iso_half; mvc_1; mvc_2; mvc_3; cc; fat];
signals = NaN(length(files));
for i = 1:1:length(signals)
    signals(i) = EMG_signal(files(i), fs);
end
%% 1)
% 1.Determine where, in the data set, the EMG signal starts and stops, so 
% that you are processing only signal that contains EMG information.

