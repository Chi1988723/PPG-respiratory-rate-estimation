function [psdResult] = myFFT(signalHandle)
% 0423 标准的FFT函数，用于再给定范围内寻找功率谱峰值
%   输入信号序列和频率范围，默认两倍补
    signal = signalHandle.signal;
    freRange = signalHandle.freRange;
    isPlot = signalHandle.isPlot;
    fs = signalHandle.fs;
    n = signalHandle.fftLenMuti;
    
    minScale = fs / n / length(signal);
    freScale = (round(freRange(1) / minScale) : round(freRange(2) / minScale));
   
    signal = signal - mean(signal);
    
    powerSpectrum = abs(fft(signal, n * length(signal)));
    [~, temp] = max(abs(powerSpectrum(round(freScale) + 1)));
    psdResult = round(freScale(temp)) * minScale;
    
    if isPlot
        figure;
%         plot(round(freScale) * minScale, powerSpectrum(round(freScale) + 1)); grid on; xlabel('频率 Hz'); title('功率谱');
        plot((0 : n * length(signal) - 1) * minScale, powerSpectrum); grid on; xlabel('频率 Hz'); title('功率谱');
%         plot(powerSpectrum); grid on; xlabel('频率 Hz'); title('功率谱');
    end
end

