close all;
clc
clear all;
% figure;


lineStyle = ['o-';  'h-';  '*-';  '.-';  'x-';  's-';  'd-';  '^-';  'p-'; '+-'; '<-'];

% signal generation;如果想要进行100组独立的测试，可以建立100次循环，产生100组独立的数据
the = zeros(1,11);
for j = 1:11  % bit per symbol: 1. BPSK; 2. QPSK; 3.8QAM; 4. 16QAM; 5. 32QAM; 6.64QAM; 7.2FSK; 8.4FSK; 9.8FSK; 10.4ASK; 11.8ASK

    if j < 10
        System.BitPerSymbol = j;
    else
        System.BitPerSymbol = j - 8;
    end
    snr = -5:30;  %SNR信噪比的设置，单位dB
    fs1 = zeros(1,16);

    acc = zeros(1, length(snr));
    for snrIndex= 1:length(snr)
        cnt = 0;
        for repeat = 1:100
            Tx.SampleRate = 32e9; %symbol Rate，信号的码元速率，可以自行定义
            Tx.Linewidth = 0;%发射信号的载波的线宽，一般与信号的相位噪声有关，大小可自行设置，这里暂时设置为0
            Tx.Carrier = 0;%发射信号的载波频率，可自行设置，这里暂设为0
            M = 2^System.BitPerSymbol;

            % 产生码元
            if j == 7 || j == 8 || j == 9 % 2fsk一个码元1位，4fsk一个码元2位
                Tx.DataSymbol = randi([0 M/64-1],1,10000);
            else
                Tx.DataSymbol = randi([0 M-1],1,10000);%每一次随机产生的数据量，这里暂时设为数据点个数为10000个
            end

            %数据的不同调制方式产生：这里把8QAM和fsk的形式单独拿出来设置，是为了实现最优的星型8QAM星座图
            if j == 7 || j == 8 || j == 9 %fsk            
%                 Tx.DataConstel = fskmod(Tx.DataSymbol,M/64,50,j,15000);
                    Tx.DataConstel = fskmod(Tx.DataSymbol, M/64, log2((M/64)), 2, 32);
            elseif j > 9 % ASK信号直接就是dataConstel，不需要另外调制
                Tx.DataConstel = Tx.DataSymbol;
            elseif M ~= 8
                h = modem.qammod('M', M, 'SymbolOrder', 'Gray');
                Tx.DataConstel = modulate(h,Tx.DataSymbol);
%                     Tx.DataConstel = qammod(Tx.DataSymbol, M, 'gray');
            else %8QAM
                    tmp = Tx.DataSymbol;
                    tmp2  = zeros(1,length(Tx.DataSymbol));
                    for kk = 1:length(Tx.DataSymbol)
                        switch tmp(kk)
                            case 0
                                tmp2(kk) = 1 + 1i;
                            case 1
                                tmp2(kk) = -1 + 1i;
                            case 2
                                tmp2(kk) = -1 - 1i;
                            case 3
                                tmp2(kk) = 1 - 1i;
                            case 4
                                tmp2(kk) = 1+sqrt(3);
                            case 5
                                tmp2(kk) = 0 + 1i .* (1+sqrt(3));
                            case 6
                                tmp2(kk) = 0 - 1i .* (1+sqrt(3));
                            case 7
                                tmp2(kk) = -1-sqrt(3);
                        end
                    end
                    Tx.DataConstel = tmp2;
                    clear tmp tmp2;
            end

            Tx.Signal = Tx.DataConstel;

            %数据的载波加载，考虑到相位噪声等
            N = length(Tx.Signal);
            dt = 1/Tx.SampleRate;
            t = dt*(0:N-1);
            Phase1 = [0, cumsum(normrnd(0,sqrt(2*pi*Tx.Linewidth/(Tx.SampleRate)), 1, N-1))];
            carrier1 = exp(1i*(2*pi*t*Tx.Carrier + Phase1));
            Tx.Signal = Tx.Signal.*carrier1;

            if(j == 1 || j >=10 ) %2PSK 和 ASK的星座图都是实数的，要转成复数
                Tx.Signal = complex(Tx.Signal);
            end

            Rx.Signal = awgn(Tx.Signal,snr(snrIndex),'measured');%数据在AWGN信道下的接收
%             Rx.Signal = Tx.Signal;

            CMAOUT = Rx.Signal;

            %normalization接收信号功率归一化
    %         CMAOUT=CMAOUT/sqrt(mean(abs(CMAOUT).^2));
    
%             subplot(1,4,snrIndex); 
%             plot(Rx.Signal,'.');

            % 检测调制类型
            type = classify(Rx.Signal);
            
            % 记录一次成功识别
            if type == j 
                cnt = cnt +1;
            else
                CMAOUT = CMAOUT/sqrt(mean(abs(CMAOUT).^2));
                type = classify2(CMAOUT, snr(snrIndex));
                if(type == j)
                    cnt = cnt + 1;
                end
            end
            
%             fs1(snrIndex) = classify(Rx.Signal);

    %         si = [real(Rx.Signal)' imag(Rx.Signal)'];
    %         center = subclust(si, [0.1 0.1]);
    %         hold on;
    %         center = complex(center( : , 1)', center(:,2)');
    %         plot(center,'*');
    %         options = statset('MaxIter',1000);
    %         gmfit = fitgmdist(si,3,'CovarianceType','diagonal',...
    %             'SharedCovariance','true','Options',options);
    %         center = cluster(gmfit, si);
    %         hold on;
    %         plot(center, '.');
        
        end
        
        % 计算该信噪比下的识别率
        acc(snrIndex) = cnt/100;

    end

%     hold on
    figure(1);
    subplot(3, 4, j);
    plot(snr, acc,lineStyle(j,:));
    axis([-5 30 0 1.1]);
    ylabel('识别正确率/%');
    xlabel('信噪比/dB');
    if(j == 1)
        title('BPSK调制方式识别');
    elseif(j == 2)
        title('QPSK调制方式识别');
    elseif(j == 3)
        title('8QAM调制方式识别');
    elseif(j == 4)
        title('16QAM调制方式识别');
    elseif(j == 5)
        title('32QAM调制方式识别');
    elseif(j == 6)
        title('64QAM调制方式识别');
    elseif(j == 7)
        title('2FSK调制方式识别');
    elseif(j == 8)
        title('4FSK调制方式识别');
    elseif(j == 9)
        title('8FSK调制方式识别');
    elseif(j == 10)
        title('4ASK调制方式识别');
    else
        title('8ASK调制方式识别');
    end
    
    
    
    figure(2);
    subplot(3, 4, j);
    plot(real(CMAOUT),imag(CMAOUT),'.'); 
    if(j == 1)
        title('BPSK调制方式星座图');
    elseif(j == 2)
        title('QPSK调制方式星座图');
    elseif(j == 3)
        title('8QAM调制方式星座图');
    elseif(j == 4)
        title('16QAM调制方式星座图');
    elseif(j == 5)
        title('32QAM调制方式星座图');
    elseif(j == 6)
        title('64QAM调制方式星座图');
    elseif(j == 7)
        title('2FSK调制方式星座图');
    elseif(j == 8)
        title('4FSK调制方式星座图');
    elseif(j == 9)
        title('8FSK调制方式星座图');
    elseif(j == 10)
        title('4ASK调制方式星座图');
    else
        title('8ASK调制方式星座图');
    end
end
grid on
legend('BPSK', 'QPSK', '8QAM', '16QAM', '32QAM', '64QAM', '2FSK', '4FSK', '8FSK', '4ASK', '8ASK');
