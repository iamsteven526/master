%3 user BPSK + LDPC

clc;
clear;


fid=fopen('H_2048_1024_z64_0635.txt','r');
C = textscan(fid,'%d','delimiter',sprintf(' '));
A = cell2mat(C);
B = reshape(A,2048,1024)';
B = sparse(double(B));

K = 1;
V = 3; % user

N = 2048; % SCMA signals in frame
M = 2; % BPSK = 2 QPSK = 4 ....
R = 0.5;  %code rate
EbN0 = 4:2:20;
SNR  = EbN0 + 10*log10(2*R*log2(M));

Nerr  = zeros(V, length(SNR));
Nbits = zeros(V, length(SNR));
BER   = zeros(V, length(SNR));

maxNumErrs = 20000;
maxNumBits = 6e8;
Niter      = 5;
ldpcDecoder = comm.LDPCDecoder(B);
ldpcEncoder = comm.LDPCEncoder(B);
ldpcDecoder.MaximumIterationCount = 250;
ldpcDecoder.IterationTerminationCondition = 'Parity check satisfied';

for k = 1:length(SNR)

    N0 = 1/(10^(SNR(k)/10)); % noise power

    while ((min(Nerr(:,k)) < maxNumErrs) && (Nbits(1,k) < maxNumBits))
        display(k);
        display(Nbits(1,k));
        %x = randi([0 M-1], V, N); % log2(M)-bit symbols
        dam = randi([0 1], V, N/2); %32400      
        for pp = 1:V
           w(pp,:) = ldpcEncoder(dam(pp,:)');%64800
           %w(pp,:) = nrPolarEncode(dam(pp,:)',N,10,false);
        end
        
        h = 1/sqrt(2)*(randn(V, N)+1j*randn(V, N)); % Rayleigh channel
        %h = 1/sqrt(2)*(repmat(randn(1, V, N), K, 1)+1j*repmat(randn(1, V, N), K, 1));
        popo = 8;
        for pp = 1:popo
            for qq = 1:N/popo
                h(:,(pp-1)*N/popo+qq) = h(:,(pp-1)*N/popo+1);
            end
        end
        for pp = 1:V
           for tt = 1:N/log2(M)
               x(pp,tt) = 1 - 2*w(pp,tt);
               %x(pp,tt) = w(pp,2*tt) + 2*w(pp,2*tt-1); %ex:10 = 2 01 = 1 %32400
           end
        end
        s = sum(h.*x,1); % joint encoding and fading channel propagation
        y = awgn(s, SNR(k));

        LLR = scmadec777(y, h, N0);
        LLR(LLR==inf) = 1500;
        LLR(LLR==-inf) = -1500;
        % symbol to bit conversion
        %r    = de2bi(x, log2(M), 'left-msb');

        for pp = 1:V
            %ansbit(pp,:) = nrPolarDecode(datar(:,pp),N/2,N,8,10,false,24);
            ansbit(pp,:) = ldpcDecoder(LLR(pp,:)');
        end

        err        = sum(xor(dam', ansbit'));
        Nerr(:,k)  = Nerr(:,k) + err.';
        Nbits(:,k) = Nbits(:,k) + log2(M)*N*R;
    end
    BER(:,k) = Nerr(:,k)./Nbits(:,k);
    k
    display(sum(Nerr(:,k))/3)
end
plot(EbN0,log10(BER))