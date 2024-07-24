% Sampling frequency
fs = 250;

% Load data
str = "new_2channels_nobiasing4.csv"; % Replace with your actual file path
T = readtable(str);             
x = T.Var2;

% Divide the signal into 4 segments of 15-second each
segment_duration = 15; % duration of each segment(second)
samples_per_segment = fs * segment_duration;
number_of_segments = floor(length(x) / samples_per_segment);

% Preallocate arrays to store SNR values and peak frequencies
SNR_values = zeros(1, number_of_segments);
peak_frequencies = zeros(1, number_of_segments);

% Preallocate array to store peak amplitudes for each segment
peak_amplitudes = zeros(1, number_of_segments);


for k = 1:number_of_segments
    % Extract the segment
    segment_start = (k-1) * samples_per_segment + 1;
    segment_end = k * samples_per_segment;
    segment = x(segment_start:segment_end);
    
    % FFT of the segment
    xdft_segment = fft(segment);
    xdft_segment = xdft_segment(1:length(segment)/2+1);
    psdx_segment = (1/(fs*length(segment)))*abs(xdft_segment).^2;
    psdx_segment(2:end-1) = 2*psdx_segment(2:end-1);
    
    % Frequency vector
    freq_segment = 0:fs/length(segment):fs/2;
    
    % Find the peak power within 12 to 13 Hz
    target_freq_indices = (freq_segment >= 12.1) & (freq_segment <= 12.9);
    [target_peak_power, target_peak_index] = max(psdx_segment(target_freq_indices));
    peak_freq = freq_segment(target_freq_indices);
    peak_freq = peak_freq(target_peak_index); % Corrected frequency
    peak_frequencies(k) = peak_freq;
    
    % Compute SNR for this segment
    % Define noise band (excluding peak +/- 0.5 Hz),
    % Data are excluded within +/- 0.5 Hz of the peak frequency
    noise_band_width = 0.5;
    noise_indices = (freq_segment >= 5) & (freq_segment <= 20) & ...
                    ~(freq_segment >= peak_freq - noise_band_width & freq_segment <= peak_freq + noise_band_width);
    noise_powers = psdx_segment(noise_indices);
    noise_power_avg = mean(noise_powers);
    SNR_values(k) = 10 * log10(target_peak_power / noise_power_avg);
    
    % Store the peak amplitude in dB/Hz
    peak_amplitudes(k) = pow2db(target_peak_power);
    
    % Plot for this segment
%     figure;
%     plot(freq_segment, pow2db(psdx_segment));
%     grid on;
%     title(sprintf("Periodogram Using FFT for Segment %d", k));
%     xlabel("Frequency (Hz)");
%     ylabel("Power/Frequency (dB/Hz)");
%     hold on;
%     plot(peak_freq, pow2db(target_peak_power), 'r*', 'MarkerSize', 10);
%     text(peak_freq, pow2db(target_peak_power), sprintf('Peak: (%0.2f Hz, %0.2f dB/Hz)', peak_freq, pow2db(target_peak_power)));
%     hold off;
end

% Display the SNR for the segements
disp('SNR for each 15-second segment:');
disp(SNR_values);

% Display the peak frequencies for all segments
disp('Peak frequency for each 15-second segment:');
disp(peak_frequencies);

% Display the peak amplitudes for all segments using disp
disp('Peak amplitudes(dB/Hz) for each 15-second segment:');
disp(peak_amplitudes);


% Perform FFT on the entire data set
xdft_full = fft(x);
xdft_full = xdft_full(1:length(x)/2+1);
psdx_full = (1/(fs*length(x)))*abs(xdft_full).^2;
psdx_full(2:end-1) = 2*psdx_full(2:end-1);

% Frequency vector for the full data set
freq_full = 0:fs/length(x):fs/2;

% Find the peak power within 12 to 13 Hz for the full data set
target_freq_indices_full = (freq_full >= 12.1) & (freq_full <= 12.9);
[target_peak_power_full, target_peak_index_full] = max(psdx_full(target_freq_indices_full));
peak_freq_full = freq_full(target_freq_indices_full);
peak_freq_full = peak_freq_full(target_peak_index_full); % Corrected frequency

% Compute SNR for the entire data set
% Define noise band (excluding peak +/- 0.5 Hz for the full data set)
noise_indices_full = (freq_full >= 5) & (freq_full <= 20) & ...
                     ~(freq_full >= peak_freq_full - 0.5 & freq_full <= peak_freq_full + 0.5);
noise_powers_full = psdx_full(noise_indices_full);
noise_power_avg_full = mean(noise_powers_full);
SNR_full = 10 * log10(target_peak_power_full / noise_power_avg_full);

% Plot the periodogram for the full data set
figure;
plot(freq_full, pow2db(psdx_full));
grid on;
title("Periodogram Using FFT for Full Data Set");
xlabel("Frequency (Hz)");
ylabel("Power/Frequency (dB/Hz)");
hold on;
plot(peak_freq_full, pow2db(target_peak_power_full), 'r*', 'MarkerSize', 10);
text(peak_freq_full, pow2db(target_peak_power_full), sprintf('Peak: (%0.2f Hz, %0.2f dB/Hz)', peak_freq_full, pow2db(target_peak_power_full)));
hold off;

% Display the SNR and peak frequency for the full data set
fprintf('Peak frequency for the entire 1-minute data: %.2f Hz\n', peak_freq_full);
fprintf('SNR for the entire 1-minute data: %.2f dB\n', SNR_full);
fprintf('Peak amplitude for the entire 1-minute data: %.2f dB/Hz\n', pow2db(target_peak_power_full));
fprintf('Average noise power (5 Hz to 20 Hz excluding peak): %.2f dB/Hz\n', pow2db(noise_power_avg_full));