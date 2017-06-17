% close all;
% clc;
figure;
% signal generation;�����Ҫ����100������Ĳ��ԣ����Խ���100��ѭ��������100�����������
the = zeros(1,9);
for j = 1:9  % bit per symbol: 1. BPSK; 2. QPSK; 3.8QAM; 4. 16QAM; 5. 32QAM; 6.64QAM;7.FSK;8:4ASK;9.8ASK
    if j < 8
        System.BitPerSymbol = j;
    else
        System.BitPerSymbol = mod(j,6);
    end
    snr = 0:100;  %SNR����ȵ����ã���λdB
%     figure;
    fs1 = zeros(1,16);
    for snrIndex= 1:length(snr)

        Tx.SampleRate = 32e9; %symbol Rate���źŵ���Ԫ���ʣ��������ж���
        Tx.Linewidth = 0;%�����źŵ��ز����߿�һ�����źŵ���λ�����йأ���С���������ã�������ʱ����Ϊ0
        Tx.Carrier = 0;%�����źŵ��ز�Ƶ�ʣ����������ã���������Ϊ0
        M = 2^System.BitPerSymbol;

        Tx.DataSymbol = randi([0 M-1],1,10000);%ÿһ�������������������������ʱ��Ϊ���ݵ����Ϊ10000��

        %���ݵĲ�ͬ���Ʒ�ʽ�����������2^3��8QAM������ʽ�����ó������ã���Ϊ��ʵ�����ŵ�����8QAM����ͼ
        if j == 7 %fsk
            Tx.DataConstel = fskmod(Tx.DataSymbol,M,50,j,15000);
        elseif j > 7            
            Tx.DataConstel = Tx.DataSymbol;            
        elseif M ~= 8
%             h = modem.qammod('M', M, 'SymbolOrder', 'Gray');
%             Tx.DataConstel = modulate(h,Tx.DataSymbol);
                Tx.DataConstel = qammod(Tx.DataSymbol, M, 'gray');
        else
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

        %���ݵ��ز����أ����ǵ���λ������
        N = length(Tx.Signal);
        dt = 1/Tx.SampleRate;
        t = dt*(0:N-1);
        Phase1 = [0, cumsum(normrnd(0,sqrt(2*pi*Tx.Linewidth/(Tx.SampleRate)), 1, N-1))];
        carrier1 = exp(1i*(2*pi*t*Tx.Carrier + Phase1));
        Tx.Signal = Tx.Signal.*carrier1;

        if(M==2)
            Tx.Signal = complex(Tx.Signal);
        end
        
%         Rx.Signal = awgn(Tx.Signal,snr(snrIndex),'measured');%������AWGN�ŵ��µĽ���
        Rx.Signal = Tx.Signal;

        CMAOUT = Rx.Signal;

        %normalization�����źŹ��ʹ�һ��
        CMAOUT=CMAOUT/sqrt(mean(abs(CMAOUT).^2));

%         subplot(1,4,snrIndex); 
%         plot(Rx.Signal,'.');
                
        fs1(snrIndex) = classify(Rx.Signal);
        
        
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
    plot(fs1,'.-');
    grid on
    hold on
    the(j) = mean(fs1);
end
% legend('QPSK', '8QAM', '16QAM', '32QAM', '64QAM','FSK','4ASK','8ASK');
legend('BPSK', 'QPSK', '8QAM', '16QAM', '32QAM', '64QAM','FSK','4ASK','8ASK');
% legend( 'QPSK', '8QAM', '16QAM', '32QAM', '64QAM');






