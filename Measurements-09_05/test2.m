 
fs = 250;
    str = "new_2channels_nobias_1.csv";
    %str = "node1_1.csv";
    T = readtable(str);             
   % t = T.Var1;
    x = T.Var2(2:end);
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
    

    dt = 1/fs; % 每个采样点之间的时间间隔

    % 创建时间向量
    t = 0:dt:(N-1)*dt;

    figure(2); 
    plot(t, x);
    
    grid on; % 在第二个图形上添加网格线
    xlabel('Time (s)');
    ylabel('EEG Signal');
    title('EEG Signal Over Time');
    
    hold on
    % 创建样本索引向量
    samples = 1:N;

%     % 绘制图形
%     figure(3);
%     plot(samples, x);
%     xlabel('Samples');
%     ylabel('EEG Signal');
%     title('EEG Signal per Sample');
%     grid on;
%     
%     
% % 设计陷波滤波器
% notch_freq = 1; % 陷波滤波器的中心频率（Hz）
% bandwidth = 0.5; % 陷波器的带宽（Hz）
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
% title('Filtered EEG Signal (1 Hz Notch Filter)');
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