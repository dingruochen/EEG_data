

function calculate_impedance
    fs = 250;  % Sampling rate in Hz
    str = 'ear_nobias1.csv';  % File name
    T = readtable(str);  % Read data from the table
    x = T.Var2; % Assuming the EEG data is in the second column

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
    
    
    target_frequency = 50; % Target frequency in Hz

    % Calculate the amplitude using Goertzel algorithm
    amplitude = goertzel(x, target_frequency, fs);

    fprintf('Amplitude at %.1f Hz: %.6f V\n', target_frequency, amplitude);

    if amplitude == 0
        fprintf('Error: Amplitude is zero. Check the input data and parameters.\n');
        return;
    end

    % Convert amplitude to Vrms assuming a sinusoidal component at this frequency
    vrms = amplitude / sqrt(2);  % Conversion for pure tone

    fprintf('Vrms at %.1f Hz: %.6f V\n', target_frequency, vrms);
    
    % Impedance calculation (adjust values as necessary for your specific setup)
    impedance = impedance_formula_ac(vrms);
    fprintf('Impedance at %.1f Hz: %.6f ohms\n', target_frequency, impedance);
end

function amplitude = goertzel(data, target_frequency, fs)
    N = length(data);
    k = round(target_frequency / fs * N);
    omega = 2 * pi * k / N;
    cosine = cos(omega);
    sine = sin(omega);
    coeff = 2 * cosine;
    q1 = 0;
    q2 = 0;
    for i = 1:N
        q0 = coeff * q1 - q2 + data(i);
        q2 = q1;
        q1 = q0;
    end
    realPart = q1 - q2 * cosine;
    imagPart = q2 * sine;
    amplitude = sqrt(realPart^2 + imagPart^2);
end

function impedance = impedance_formula_ac(vrms)
    % Example constants, adjust according to actual application
    gain = 24.0; % Amplifier gain used for measurement
    current_amplitude_nano_amp = 6.0; % Current amplitude in nanoamps
    irms = current_amplitude_nano_amp * 1e-9 / sqrt(2); % Irms calculation
    impedance = abs(vrms / irms - 6800.0); % Example internal impedance calculation
end


    



% fs = 250;
%     str = "ear_nobias2.csv";
%     %str = "ROPT_10K_1_21_3"
%     T = readtable(str);             
%    % t = T.Var1;
%     x = T.Var2;
%     N = length(x);
%     xdft = fft(x);
%     xdft = xdft(1:N/2+1);
%     psdx = (1/(fs*N))*abs(xdft).^2;
%     psdx(2:end-1) = 2*psdx(2:end-1);
%     freq = 0:fs/length(x):fs/2;
%     figure(1)
%     plot(freq,pow2db(psdx))
%     
%     grid on
%     title("Periodogram Using FFT")
%     xlabel("Frequency (Hz)")
%     ylabel("Power/Frequency (dB/Hz)")
%     hold on
%     P_I = pow2db(psdx)
%     
% 
% %     dt = 1/fs; % 每个采样点之间的时间间隔
% % 
% %     % 创建时间向量
% %     t = 0:dt:(N-1)*dt;
% 
% %     figure(2); 
% %     plot(t, x);
% %     xlabel('Time (s)');
% %     ylabel('EEG Signal');
% %     title('EEG Signal Over Time');
% %     grid on; % 在第二个图形上添加网格线
%     
%     % 创建样本索引向量
%     samples = 1:N;
% 
%     % 绘制图形
%     figure(3);
%     plot(samples, x);
%     xlabel('Samples');
%     ylabel('EEG Signal');
%     title('EEG Signal per Sample');
%     grid on;
    
    
% % 设计陷波滤波器
% notch_freq = 50; % 陷波滤波器的中心频率（Hz）
% bandwidth = 10; % 陷波器的带宽（Hz）
% wo = notch_freq/(fs/2); % 归一化频率
% bw = bandwidth/(fs/2); % 归一化带宽
% [b, a] = iirnotch(wo, bw);
% 
% % 应用滤波器
% filtered_x = filter(b, a, x);
% 
% 
% 
% % 绘制滤波后的信号
% figure(4); % 新建一个图形窗口
% plot(samples, filtered_x);
% xlabel('Samples');
% ylabel('Filtered EEG Signal');
% title('Filtered EEG Signal (48-52 Hz Notch Filter)');
% grid on;
% 
% 
% 
% % 对滤波后的信号执行FFT
% filtered_xdft = fft(filtered_x);
% filtered_xdft = filtered_xdft(1:N/2+1);
% 
% % 计算滤波后信号的PSD
% filtered_psdx = (1/(fs*N)) * abs(filtered_xdft).^2;
% filtered_psdx(2:end-1) = 2*filtered_psdx(2:end-1);
% 
% % 绘制滤波后信号的频率响应
% figure(5); % 新建一个图形窗口
% plot(freq, pow2db(filtered_psdx));
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
% title('Periodogram of Filtered Signal');
% grid on;