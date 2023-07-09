function [fullAcfResult] = myHarmonicPower(signalHandle)
% 0624 

    signal = signalHandle.signal;
    fs = signalHandle.fs;
    HRfs = signalHandle.HRfs;
    periodRange = signalHandle.periodRange;
    freRange = flipud(1 ./ periodRange);
    hbRange = signalHandle.HRfreRange;
    HRsignal = signalHandle.HRsignal;
    %% Power Spectrum
    signal = signal - mean(signal);
    nFFT = length(signal) * 20;
    powerSpectrum = fft(signal .* hamming(length(signal)) , nFFT);%
    powerSpectrum = powerSpectrum .* conj(powerSpectrum);
    
    HRsp = abs(fft(HRsignal));
    nFFThr = length(HRsp);
    hbScale = round(hbRange(1) / HRfs * nFFThr : hbRange(2) / HRfs * nFFThr);
    [~, temp] = max(HRsp(hbScale));
    hbFre = hbScale(temp) / nFFThr * HRfs * 0.5;
    freRange(2) = min(hbFre, freRange(2));
    
    powerSpectrum = powerSpectrum - mean(powerSpectrum(floor(freRange(1) / fs * nFFT) : floor(hbFre / fs * nFFT)));
    %% 和谱，直接功率谱求平均计算
    delayScale = (floor(freRange(1) / fs * nFFT) : floor(freRange(2) / fs * nFFT));
    sumSpectrum = zeros(length(delayScale), 1);
    for ii = 1 : length(delayScale)
        sumSpectrum(ii) = sum(powerSpectrum(delayScale(ii) + 1 : delayScale(ii) : (floor(hbFre / fs * nFFT)))); % 求和范围由心跳决定
    end
       
    %% 0625 倍数周期判断
    [~, temp] = max(sumSpectrum); 
    if isempty(temp)
        fullAcfResult = NaN;
    else
        fullAcfResult = delayScale(temp) * fs / nFFT;
    end
    
    %% 决定是否画图
    if signalHandle.isPlot
        figure; hold on;
        plot( (1 * delayScale * fs ./ length(powerSpectrum)), sumSpectrum,  'linewidth', 1);xlabel('Respiratory rate','fontsize', 12);ylabel('Amplitude','fontsize', 12);
    end
end

