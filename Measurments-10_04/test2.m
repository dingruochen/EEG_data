 
    fs = 250;
    str = "COMP_EAR1_10_4.csv";
    T = readtable(str);             
    x = T.Var2;
    N = length(x);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(fs*N))*abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:fs/length(x):fs/2;
    figure(1)
    plot(freq,pow2db(psdx))
    
    grid on
    title("Periodogram Using FFT")
    xlabel("Frequency (Hz)")
    ylabel("Power/Frequency (dB/Hz)")
    hold on
    P_I = pow2db(psdx)
    

    samples = 1:N;

    figure(3);
    plot(samples, x);
    xlabel('Samples');
    ylabel('EEG Signal');
    title('EEG Signal per Sample');
    grid on;
    
 
% notch_freq = 50;
% bandwidth = 10; 
% wo = notch_freq/(fs/2);
% bw = bandwidth/(fs/2); 
% [b, a] = iirnotch(wo, bw);
% 
% filtered_x = filter(b, a, x);
% 
% 
% 
%
% figure(4); 
% plot(samples, filtered_x);
% xlabel('Samples');
% ylabel('Filtered EEG Signal');
% title('Filtered EEG Signal (48-52 Hz Notch Filter)');
% grid on;
% 
% 
% 
% % FFT
% filtered_xdft = fft(filtered_x);
% filtered_xdft = filtered_xdft(1:N/2+1);
% 
% % PSD
% filtered_psdx = (1/(fs*N)) * abs(filtered_xdft).^2;
% filtered_psdx(2:end-1) = 2*filtered_psdx(2:end-1);
% 
%
% figure(5);
% plot(freq, pow2db(filtered_psdx));
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
% title('Periodogram of Filtered Signal');
% grid on;