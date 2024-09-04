% Sample file list, replace with your actual filenames
fileList = ["bias_usb3.csv"]; 

fs = 250; % Sampling frequency

% Loop over each file in the list
for i = 1:length(fileList)
    % Read the table from the file
    T = readtable(fileList(i), 'VariableNamingRule', 'preserve');
    x = T{:, 2}; % Access the second column
    N = length(x);
    
    % Compute the FFT
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    
    % Compute the Power Spectral Density (PSD)
    psdx = (1/(fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    
    % Frequency vector
    freq = 0:fs/length(x):fs/2;
    
    % Convert PSD to decibels
    P_I = pow2db(psdx);
    
    % Create a new figure for each file
    figure;
    
    % Plot the PSD
    plot(freq, P_I);
    
    % Customize the plot
    grid on;
    title(sprintf("Periodogram Using FFT for File %d: %s", i, fileList(i)));
    xlabel("Frequency (Hz)");
    ylabel("Power/Frequency (dB/Hz)");
end
