
fs = 250;
    str = "2elec_ear_far_ssvep3.csv";
    T = readtable(str);             
   % t = T.Var1;
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
    
    
