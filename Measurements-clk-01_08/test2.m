fs = 250;
str = "2channels_new3.csv";
T = readtable(str, 'VariableNamingRule', 'preserve'); % Preserve original column names
x = T{:, 3}; % Access the second column
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(fs*N))*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:fs/length(x):fs/2;
figure(1)
plot(freq, pow2db(psdx))

grid on
title("Periodogram Using FFT")
xlabel("Frequency (Hz)")
ylabel("Power/Frequency (dB/Hz)")
hold on
P_I = pow2db(psdx)