% form 2D matrices into 3D matrices
st3d = reshape(sig2d, num_rows, num_columns, 512);
sig3d = reshape(sig2d, num_rows, num_columns, 512);

% calculate averaged spectra for every 10 second section (every row)
% apply Fast Fourier Transform
spectrums = [];
for i = range(num_rows)
    sig3d = sig3d(i,:,:);
    [f_axis, fourier_t] = sp.fast_fourier_t(sig3d, fs);
    % calculate Fourier power spectrum
    fourier_t = mean(fourier_t.^2);
    spectrums(end+1) = fourier_t;
end