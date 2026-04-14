function [timeVec, freqVec, tf_power] = multitaper_tfr(data, fs, freq_range, window_length, step_size, time_bw)
% Multitaper time-frequency power (like FieldTrip's mtmconvol)
%
% Inputs:
%   data         - 1D signal (1 x time)
%   fs           - sampling rate (Hz)
%   freq_range   - [fmin fmax] frequency range (Hz)
%   window_length- analysis window length (s)
%   step_size    - time step between windows (s)
%   time_bw      - time-bandwidth product (e.g., 3)
%
% Outputs:
%   timeVec      - time centers of each window
%   freqVec      - frequency vector
%   tf_power     - [freq x time] power matrix

win_samples = round(window_length * fs);
step_samples = round(step_size * fs);
nfft = 2^nextpow2(win_samples);

% Define frequency vector
f = fs * (0:(nfft/2)) / nfft;
freq_idx = find(f >= freq_range(1) & f <= freq_range(2));
freqVec = f(freq_idx);

% Multitaper parameters
K = floor(2 * time_bw - 1);
[tapers, ~] = dpss(win_samples, time_bw, K);

% Slide window across signal
nSteps = floor((length(data) - win_samples) / step_samples) + 1;
tf_power = zeros(length(freq_idx), nSteps);
timeVec = zeros(1, nSteps);

for i = 1:nSteps
    start_idx = (i - 1) * step_samples + 1;
    seg = data(start_idx : start_idx + win_samples - 1);
    timeVec(i) = (start_idx + win_samples/2) / fs;

    % Taper & FFT
    Xk = zeros(K, nfft);
    for k = 1:K
        tapered = seg .* tapers(:, k)';
        Xk(k, :) = fft(tapered, nfft);
    end

    % Average power
    Pk = abs(Xk(:, 1:nfft/2+1)).^2;
    power_avg = mean(Pk, 1) / fs;

    tf_power(:, i) = power_avg(freq_idx)';
end
end