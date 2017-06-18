function [type] = classify(x)
    type = -1;
%     1. BPSK; 2. QPSK; 3.8QAM; 4. 16QAM; 5. 32QAM; 6.64QAM; 7.2FSK; 
%     8.4FSK; 9.8FSK;  10.4ASK; 11.8ASK
%     type1th
%     272.0000   33.9599   43.5059   13.9654    1.9869   11.5302    0.5976
%     0.0469    0.0106  48.2447   66.1035

    th1 = [272.0000   33.9622   43.5034   13.9617 11.5009  48.2853   66.0788 0];
%     1. BPSK; 2. QPSK; 3.8QAM; 4. 16QAM; 5.64QAM 6.4ASK; 7.8ASK 8.MFSK&32QAM

    th2 = [15.8673    0.8983    3.7164]; %1.2FSK 2.4FSK 3.8FSK 4.32QAM 6.7745

    th3 = [3.2305 0];%1.32QAM 2.8FSK

    m80 = MMk(x, 8, 0);
    m60 = MMk(x, 6, 0);
    m40 = MMk(x, 4, 0);
    m20 = MMk(x, 2, 0);
    m21 = MMk(x, 2, 1);
    
    c80 = m80  - 28*m20*m60 - 35*(m40^2) + 420*(m20^2)*m40 - 630*(m20^4);
    c21 = m21;
    
    fs1 = abs(c80)/(abs(c21)^4);
    
    diff1 = abs(th1-fs1);
    fil1 = find(diff1 == min(diff1));
    
    if fil1 ~= 8
        switch fil1
            case 1
                type = 1;
            case 2
                type = 2;
            case 3
                type = 3;
            case 4
                type = 4;
            case 5
                type = 6;
            case 6
                type = 10;
            case 7
                type = 11;
        end
    else
        xt = real(x);
        yt = imag(x);
        instantAmplitude = sqrt(xt.^2 + yt.^2);
        AmplitudeMean = mean( ( instantAmplitude - 1 ) );  
        diff3 = abs(th3 - AmplitudeMean);
        fil3 = find(diff3 == min(diff3));
        if fil3 == 1
            type = 5;
        else
            xdm = diff(x);

            md42 = MMk(xdm, 4, 2);
            md20 = MMk(xdm ,2, 0);
            md21 = MMk(xdm, 2, 1);
            md63 = MMk(xdm, 6, 3);
            md41 = MMk(xdm, 4, 1);

            cd42 = md42 - abs(md20)^2 - 2*md21^2;
            cd63 = md63 - 6*md41*md20 - 9*md42*md21 + 18*md20^2*md21 + 12*md21^3;

            fs2 = (abs(cd63)^2)/(abs(cd42)^3);

            diff2 = abs(th2-fs2);
            fil2 = find(diff2 == min(diff2));
            
            type = fil2 + 6;
            
        end        
    end
end