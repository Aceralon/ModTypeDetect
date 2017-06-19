function type = classify2(x, snr)
% 通过聚类来识别PSK, QAM
    type = -1;
%     1. BPSK; 2. QPSK; 3.8QAM; 4. 16QAM; 5. 32QAM; 6.64QAM; 

    if snr < 10
        range = 0.5;
    else
        range = 0.2;
    end


    xn = [real(x') imag(x')];
    C = subclust(xn, range);
    
    num = length(C);
    
    if(num == 2)
        type = 1;
    elseif(num == 4)
        type = 2;
    elseif(num == 8)
        type = 3;
    elseif(num == 16)
        type = 4;
    end
    
  
end