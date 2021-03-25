% Saint Petersburg Electrotechnical University, Saint Petersburg, Russia
% Faculty of Radio Engineering
% Department of Theoretical Fundamentals of Radio Engineering
% Vyacheslav P. Klimentyev and Alexander B. Sergienko, 2015

% Codebooks
clc;
clear;
CB(:,:,1) = [ 0                  0                  0                  0;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071 ];

CB(:,:,2) = [ 0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0 ];

CB(:,:,3) = [0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
              0                  0                  0                  0 ];

CB(:,:,4) = [ 0                  0                  0                  0;...
              0                  0                  0                  0;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071 ];

CB(:,:,5) = [0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
              0                  0                  0                  0;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071 ];

CB(:,:,6) = [ 0                  0                  0                  0;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0 ];
          

CB = 0.8163*CB;

K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)

N = 500; % SCMA signals in frame

EbN0 = 15:2.5:15;
SNR  = EbN0 + 10*log10(log2(M)*V/K);

Nerr  = zeros(V, length(SNR));
Nbits = zeros(V, length(SNR));
BER   = zeros(V, length(SNR));

maxNumErrs = 60000;
maxNumBits = 5e8;
Niter      = 6;

for k = 1:length(SNR)

    N0 = 1/(10^(SNR(k)/10)); % noise power

    while ((min(Nerr(:,k)) < maxNumErrs) && (Nbits(1,k) < maxNumBits))

        x = randi([0 M-1], V, N); % log2(M)-bit symbols

        h = 1/sqrt(2)*(randn(K, V, N)+1j*randn(K, V, N)); % Rayleigh channel
        %h = 1/sqrt(2)*(repmat(randn(1, V, N), K, 1)+1j*repmat(randn(1, V, N), K, 1));
        for qq = 1:N
            h(:,:,qq) = h(:,:,1);
        end
        s = scmaenc(x, CB, h); % joint encoding and fading channel propagation
        y = awgn(s, SNR(k));

        LLR = scmadec(y, CB, h, N0, Niter);

        % symbol to bit conversion
        r    = de2bi(x, log2(M), 'left-msb');
        data = zeros(log2(M)*N, V);
        for kk = 1:V
            data(:,kk) = reshape(downsample(r, V, kk-1).',[],1);
        end

        % LLR to bit conversion
        datadec = reshape((LLR <= 0), [log2(M) N*V]).';
        datar   = zeros(log2(M)*N, V);
        for kk = 1:V
            datar(:,kk) = reshape(downsample(datadec, V, kk-1).', [], 1);
        end

        err        = sum(xor(data, datar));
        Nerr(:,k)  = Nerr(:,k) + err.';
        Nbits(:,k) = Nbits(:,k) + log2(M)*N
    end
    BER(:,k) = Nerr(:,k)./Nbits(:,k);
    k
end
plot(EbN0,log10(BER))