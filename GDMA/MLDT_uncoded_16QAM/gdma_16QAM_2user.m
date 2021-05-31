%3 user BPSK + LDPC

clc;
clear;

K = 1;
V = 2; % user

N = 20; % SCMA signals in frame
M = 16; % BPSK = 2 QPSK = 4 ....
R = 1;  %code rate
EbN0 = 0:5:40;
SNR  = EbN0 + 10*log10(R*log2(M));

Nerr  = zeros(V, length(SNR));
Nbits = zeros(V, length(SNR));
BER   = zeros(V, length(SNR));

maxNumErrs = 10000;
maxNumBits = 6e8;
Niter      = 5;

for k = 1:length(SNR)

    N0 = 1/(10^(SNR(k)/10)); % noise power

    while ((min(Nerr(:,k)) < maxNumErrs) && (Nbits(1,k) < maxNumBits))
        display(k);
        display(Nbits(1,k));
        %x = randi([0 M-1], V, N); % log2(M)-bit symbols
        w = randi([0 1], V, N); %32400      

        
        h = 1/sqrt(2)*(randn(V, N/4)+1j*randn(V, N/4)); % Rayleigh channel
        popo = 1;
        for pp = 1:popo
            for qq = 1:N/(popo*4)
                h(:,(pp-1)*N/(popo*4)+qq) = h(:,(pp-1)*N/(popo*4)+1);
            end
        end

        for pp = 1:V
           x(pp,:) = qammod(w(pp,:)',16,'UnitAveragePower',true,'InputType','bit')';
        end
        s = sum(h.*x,1); % joint encoding and fading channel propagation
        y = awgn(s, SNR(k));

        LLR = scmadec777_2user(y, h, N0);
        LLR(LLR==inf) = 1500;
        LLR(LLR==-inf) = -1500;       
        LLR(LLR>= 0) = 0;
        LLR(LLR< 0) = 1;
        
        % symbol to bit conversion
        %r    = de2bi(x, log2(M), 'left-msb');

        err        = sum(xor(w', LLR'));
        Nerr(:,k)  = Nerr(:,k) + err.';
        Nbits(:,k) = Nbits(:,k) + N;
    end
    BER(:,k) = Nerr(:,k)./Nbits(:,k);
    k
    display(sum(Nerr(:,k))/2)
end
plot(EbN0,log10(BER))