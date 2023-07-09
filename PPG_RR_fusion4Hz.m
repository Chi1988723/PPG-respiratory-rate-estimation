

% 0624 提交版本


clear;
format compact;
close all;

load('CapnoData');

nData = 42;
MAE = zeros(nData, 11);
windowRetentionRatio = zeros(11, 1);
deSample = 1;
deSample2 = 75;
fsDeSampled = fs / deSample;
fsDeSampled2 = fsDeSampled / deSample2;
lenSection2 = 32 * fsDeSampled2;
lenOverlap2 = 29 * fsDeSampled2;
lenSection = 32 * fs;
lenOverlap = 29 * fs;
signalHandle.HRfreRange = [0.5; 2.5];
signalHandle.isPlot = 0;


% lpFilt1 = designfilt('lowpassfir', 'PassbandFrequency', 1.5 / fs * 2, ...% 0.4 * 
%                   'StopbandFrequency', 2 / fs * 2, 'PassbandRipple', 3, ...
%                   'StopbandAttenuation', 30);
lpFilt2 = designfilt('lowpassfir', 'PassbandFrequency', 1.1 / fsDeSampled * 2, ...% 0.4 * 
                  'StopbandFrequency', 1.5 / fsDeSampled * 2, 'PassbandRipple', 1, ...
                  'StopbandAttenuation', 30);
countThe = (-1 : 9)' ; 
nWindow = 0;
tic
for ii = 1 : 42
    fileName = manifest{(ii - 1) * 8 + 1}; % 读取mat的名字;
    load(['C:\Resources\202206\PPGRR\dataverse_files\', fileName(1 : 22)]);
    %% karlen 方法的标准MSE， 与co2范围内的平均值对比，16秒窗口，8秒一更新
    rrx = reference.rr.co2.x;
    rry = reference.rr.co2.y;
    rry(isinf(rry)) = 0;
    rrsfx = SFresults.x;
    pl = signal.pleth.y;    
    
    rrySmooth = zeros(size(rrsfx));
    for jj = 1 : 150
        timeRange = [rrsfx(jj) - 16; rrsfx(jj) + 16];
        rrySmooth(jj) = sum(rry((rrx > timeRange(1)) & (rrx < timeRange(2)))) / sum((rrx > timeRange(1)) .* (rrx < timeRange(2)));
    end
    %% 低通滤波减轻心跳的干扰：估计每个信号的心跳频率，并对应地将其滤除 
%     pl = filter(lpFilt1, pl);
%     pl = pl(1 : deSample : end);
    plLowpass = filter(lpFilt2, pl);
    plLowpass = plLowpass(1 : deSample2 : end);
    
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
%     temp = abs(rrySmooth' - amdfFre( : , 1));
    temp = abs(rrySmooth(1 : 149)' - amdfFre(2 : 150, 1));
    nWindow = nWindow + sum(~isnan(rrySmooth(1 : 149)));
%     temp = [nan; nan; abs(rrySmooth(1 : 148)' - amdfFre(3 : 150, 1))];
%     temp = [nan; nan; nan; abs(rrySmooth(1 : 147)' - amdfFre(4 : 150, 1))];
    for cc = 1 : 11
        MAE(ii, cc) = (mean(temp(~isnan(temp) & (countJudge(2 : 150) > countThe(cc))) .^ 1));
        windowRetentionRatio(cc) = windowRetentionRatio(cc) + sum((~isnan(rrySmooth(1 : 149)')) & (countJudge(2 : 150) > countThe(cc)));
    end   
    
    
%     temp = rrsfx(1 : 149);
%     temp2 = amdfFre(2 : 150, 1);
%     the = 7;
%     figure(2); clf; colororder({'k', '#EDB120'}); yyaxis left;box off; 
%     plot(rrsfx, rrySmooth,'--', 'linewidth', 2, 'color', 'black'); hold on; %title([num2str(ii),'平滑参考呼吸频率']);
% %     plot(temp, temp2,  'linewidth', 1, 'color', [0 0.4470 0.7410]);
%     figure(2); hold on; scatter(temp((countJudge(2 : 150) > countThe(the))), temp2((countJudge(2 : 150) > countThe(the))),'marker', '*', 'linewidth', 1, 'markeredgecolor', [0 0.4470 0.7410]);
%     scatter(temp((countJudge(2 : 150) < countThe(the + 1))), temp2((countJudge(2 : 150) < countThe(the + 1))),'marker', '*', 'linewidth', 1, 'markeredgecolor', [0.8500 0.3250 0.0980]);
%     xlabel('Time (second)', 'fontsize', 12); ylabel('breaths per minutes', 'fontsize', 12);
%     yyaxis right; ylim([0; 15]);
%     plot(temp, countJudge(2 : 150), 'linewidth', 1, 'color', '#EDB120', 'linestyle', ':');
%     hold off; grid off;  legend('Reference', 'Retained windows', 'Discarded windows','Counter value','fontsize', 12);box off;


end
toc
figure; hold on;
MAE = sort(MAE);
yyaxis left; ylim([0; 1.2]); ylabel('MAE (breath per minute)', 'fontsize', 12); xlabel('Counting threshold','fontsize', 12);
plot(countThe + 1, MAE(11, : ), 'linewidth', 1, 'color',[0 0.4470 0.7410], 'marker', '+');
plot(countThe + 1, mean(MAE(20 : 21, : )), 'linewidth', 1, 'color',[0 0.4470 0.7410], 'marker', 'o');
plot(countThe + 1, MAE(32, : ), 'linewidth', 1, 'color',[0 0.4470 0.7410], 'marker', '*');
yyaxis right; ylim([0; 1.2]); y = ylabel({'Window retention ratio'}, 'fontsize', 12); set(y, 'Rotation', 270, 'Position',[10.9945161290323,0.550000536441814,-1]);
plot(countThe + 1, windowRetentionRatio / nWindow, 'linewidth', 1, 'marker', 's');
legend('MAE 25th percentiles','MAE median','MAE 75th percentiles', 'Window retention ratio','fontsize', 12);


