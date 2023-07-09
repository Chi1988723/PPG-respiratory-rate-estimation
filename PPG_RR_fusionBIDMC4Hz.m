


clear;
format compact;
close all;

load('bidmc_data.mat');
fs = 125;

nData = 53;
MAE = zeros(nData, 11);
windowRetentionRatio = zeros(11, 1);
nWindow =0;
deSample = 1;
deSample2 = 31.25;
fsDeSampled = fs / deSample;
fsDeSampled2 = fsDeSampled / deSample2;
lenSection2 = 32 * fsDeSampled2;
lenOverlap2 = 29 * fsDeSampled2;
lenSection = 32 * fsDeSampled;
lenOverlap = 29 * fsDeSampled;
signalHandle.HRfreRange = [0.5; 2.5];
signalHandle.isPlot = 0;


% 0827 统一设计峰值检测中PAMDF的低通滤波器
lpFilt2 = designfilt('lowpassfir', 'PassbandFrequency', 1.1 / fsDeSampled * 2, ...% 0.4 * 
                  'StopbandFrequency', 1.5 / fsDeSampled * 2, 'PassbandRipple', 1, ...
                  'StopbandAttenuation', 30);
countThe = (-1 : 9)' ;          

nWindow = 0;
tic
for ii = 1 : nData
    dataii = data(ii);
    pl = dataii.ppg.v;
    pl = pl - mean(pl);
    ref1 = sort(dataii.ref.breaths.ann1 / fs);
    ref2 = sort(dataii.ref.breaths.ann2 / fs);

    % 人工标注的数据内有标记重复的，这可能导致在求呼吸速率时出现Inf
    ref1(diff(ref1) == 0) = [];
    ref2(diff(ref2) == 0) = [];
    
    
    ref1Time = ref1(1 : end - 1);
    ref2Time = ref2(1 : end - 1);
    ref1RR = 60 ./ diff(ref1);
    ref2RR = 60 ./ diff(ref2);
    
    rrySmooth = zeros(150, 1); % 仍然是480s的数据，采用32s滑窗，3s的覆盖，那么共150个估计窗口
    timeScale = 16 : 3 : 463;
    for jj = 1 : 150
        timeRange = [timeScale(jj) - 16; timeScale(jj) + 16];
        ref1Mean = sum(ref1RR .* (ref1Time > timeRange(1)) .* (ref1Time < timeRange(2))) / sum((ref1Time > timeRange(1)) .* (ref1Time < timeRange(2)));
        ref2Mean = sum(ref2RR .* (ref2Time > timeRange(1)) .* (ref2Time < timeRange(2))) / sum((ref2Time > timeRange(1)) .* (ref2Time < timeRange(2)));
        if abs(ref1Mean - ref2Mean) < 2
            rrySmooth(jj) = (ref1Mean + ref2Mean) / 2;
        else
            rrySmooth(jj) = nan;
        end
    end
    %% 低通滤波减轻心跳的干扰：估计每个信号的心跳频率，并对应地将其滤除 !!! 滤波器函数可能有问题 ！！！
    plLowpass = filter(lpFilt2, pl);
    plLowpass = plLowpass(round(1 : deSample2 : length(plLowpass)));
    
    %% 
    nSection = floor((length(plLowpass) - lenSection2) / (lenSection2 - lenOverlap2) + 1);
    amdfFre = zeros(1 + nSection, 3);
    signalHandle.fs = fsDeSampled2;    
    signalHandle.HRfs = fsDeSampled;    
    P = 999;
    R = 5; 
    H = 1;
    Q = 1;
    defaultPeriod = [0.01; 15];
    signalHandle.periodRange = defaultPeriod;
    
    count = 0; % 序贯估计计数
    countJudge = zeros(nSection, 1);
    for mm = 1 : nSection
        signalHandle.signal = plLowpass(round((mm - 1) * (lenSection2 - lenOverlap2) + 1 : (mm - 1) * (lenSection2 - lenOverlap2) + lenSection2));
        signalHandle.HRsignal = pl(round((mm - 1) * (lenSection - lenOverlap) + 1 : (mm - 1) * (lenSection - lenOverlap) + lenSection));
        temp = 60 * myHarmonicPower(signalHandle);  
        xt = amdfFre(mm, 1);

        signalHandle.periodRange = defaultPeriod;
        temp2 = 60 * myHarmonicPower(signalHandle);
        
        if (abs(temp - temp2) > 1) % 如果相比大范围的估计，也足够准确，那么可以
            count = count - 1;
        else
            count = min(10, count + 1);
        end
        countJudge(mm) = count; % 纪录每一次的质量

        Ptemp = P + Q;
        
        K = Ptemp * H' / (H * Ptemp * H' + R);
        P = (1 - K * H) * Ptemp;         
        amdfFre(mm + 1, 1) = xt + K * (temp - H * xt);
        div = 1 * (P / R);
        if (isnan(amdfFre(mm + 1, 1)) )
            signalHandle.periodRange = defaultPeriod;
            P = Ptemp;
            amdfFre(mm + 1, 1) = amdfFre(mm - 1, 1);
        elseif ( count < 1)
            signalHandle.periodRange = defaultPeriod;
            P = 1000;
        else
            signalHandle.periodRange = [max(60 / amdfFre(mm + 1, 1) / (1 + div), defaultPeriod(1)); min(60 / amdfFre(mm + 1, 1) * (1 + div), defaultPeriod(2))];
        end

    end
    amdfFre = amdfFre(2 : end, : );
    
    temp = abs(rrySmooth(1 : 149) - amdfFre(2 : 150, 1));
    nWindow = nWindow + sum(~isnan(rrySmooth(1 : 149)));
    for cc = 1 : 11
        MAE(ii, cc) = (mean(temp(~isnan(temp) & (countJudge(1 : 149) > countThe(cc))) .^ 1));
        windowRetentionRatio(cc) = windowRetentionRatio(cc) + sum((~isnan(rrySmooth(1 : 149))) & (countJudge(2 : 150) > countThe(cc)));
    end
    
    
%     figure(2); plot(timeScale, rrySmooth,'--', 'linewidth', 2, 'color', 'black'); hold on; grid on; %title([num2str(ii),'平滑参考呼吸频率']); 
%     figure(2); hold on; plot(timeScale, amdfFre( : , 1),'-*', 'linewidth', 1, 'color', [0 0.4470 0.7410]); grid on; xlabel('Time (second)'); ylabel('breaths / min');% legend('Sum AMDF', 'reference RR','SmartFusion', 'Fusion');
%     hold off; grid off; box off; legend('Reference', 'BW', 'AM and FM', 'Modified Kalman Fusion');

%     temp = timeScale(1 : 149);
%     temp2 = amdfFre(2 : 150, 1);
%     the = 7;
%     figure(2); clf; plot(timeScale, rrySmooth,'--', 'linewidth', 2, 'color', 'black'); hold on; %title([num2str(ii),'平滑参考呼吸频率']);
% %     plot(temp, temp2,  'linewidth', 1, 'color', [0 0.4470 0.7410]);
%     figure(2); hold on; scatter(temp((countJudge(2 : 150) > countThe(the))), temp2((countJudge(2 : 150) > countThe(the))),'marker', '*', 'linewidth', 1, 'markeredgecolor', [0 0.4470 0.7410]);
%     scatter(temp((countJudge(2 : 150) < countThe(the + 1))), temp2((countJudge(2 : 150) < countThe(the + 1))),'marker', '*', 'linewidth', 1, 'markeredgecolor', [0.8500 0.3250 0.0980]);
%     grid on; xlabel('Time (second)', 'fontsize', 12); ylabel('breaths per minutes', 'fontsize', 12);
%     
%     hold off; grid off; box off; legend('Reference', 'Retained windows', 'Discarded windows','fontsize', 12);
end
toc

figure; hold on;
MAE = sort(MAE);
yyaxis left; ylabel('MAE (breath per minute)', 'fontsize', 12); xlabel('Counting threshold','fontsize', 12);
plot(countThe + 1, (MAE(13, : ) +3 *  MAE(14, : )) / 4, 'linewidth', 1, 'color',[0 0.4470 0.7410], 'marker', '+');
plot(countThe + 1,(MAE(27, : )), 'linewidth', 1, 'color',[0 0.4470 0.7410], 'marker', 'o');
plot(countThe + 1, (MAE(40, : ) + 3 * MAE(41, : )) / 4, 'linewidth', 1, 'color',[0 0.4470 0.7410], 'marker', '*');
yyaxis right; ylim([0; 1.2]); y = ylabel({'Window retention ratio'}, 'fontsize', 12); set(y, 'Rotation', 270, 'Position',[10.9945161290323,0.50000536441814,-1]);
plot(countThe + 1, windowRetentionRatio / nWindow, 'linewidth', 1, 'marker', 's');
legend('MAE 25th percentiles','MAE median','MAE 75th percentiles', 'Window retention ratio','fontsize', 12);

