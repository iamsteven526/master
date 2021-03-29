% Saint Petersburg Electrotechnical University, Saint Petersburg, Russia
% Faculty of Radio Engineering
% Department of Theoretical Fundamentals of Radio Engineering
% Vyacheslav P. Klimentyev and Alexander B. Sergienko, 2015
% 00 01 10 11
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

          
CB = CB*0.8163;

          
K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)

N = 2048; % SCMA signals in frame
R = 0.5;
EbN0 = 10:2:10;
SNR  = EbN0 + 10*log10(R*log2(M)*V/K);   %noise power maybe wrong!!!

Nerr  = zeros(V, length(SNR));
NBLerr  = zeros(length(SNR));
Nbits = zeros(V, length(SNR));
Nblocks = zeros(length(SNR));
BER   = zeros(V, length(SNR));
BLER   = zeros(length(SNR));

maxNumErrs = 3000;
maxNumBits = 1e6;
Niter      = 8;

fid=fopen('H_2048_1024_z64_0635.txt','r');
C = textscan(fid,'%d','delimiter',sprintf(' '));
A = cell2mat(C);
B = reshape(A,2048,1024)';
B = sparse(double(B));

ldpcDecoder = comm.LDPCDecoder(B);
ldpcDecoder.MaximumIterationCount = 250;
ldpcDecoder.IterationTerminationCondition = 'Parity check satisfied';
ldpcEncoder = comm.LDPCEncoder(B);
for k = 1:length(SNR)

    N0 = 1/(10^(SNR(k)/10)); % noise power

    while ((min(Nerr(:,k)) < maxNumErrs) && (Nbits(1,k) < maxNumBits))

        %x = randi([0 M-1], V, N/2); % log2(M)-bit symbols
        display(k);
        display(Nbits(1,k))
        dam = randi([0 1], V, N/2); %32400
        for pp = 1:V
           w(pp,:) = ldpcEncoder(dam(pp,:)');%64800
           %w(pp,:) = nrPolarEncode(dam(pp,:)',N,10,false);
        end
        h = 1/sqrt(2)*(randn(K, V, N)+1j*randn(K, V, N)); % Rayleigh channel
        multiblock = 1;
        for pp = 1:multiblock
            for qq = 1:N/multiblock
                h(:,:,(pp-1)*N/multiblock+qq) = h(:,:,(pp-1)*N/multiblock+1);
            end
        end
        
        for kkk = 1:K
           for vvv = 1:V
               TimeVaryingCIR = gjmi();
               h(kkk,vvv,:) = h(kkk,vvv,:).*reshape(abs(TimeVaryingCIR),1,1,N)*(9/8);
           end
        end
        
        %h = 1/sqrt(2)*(repmat(randn(1, V, N), K, 1)+1j*repmat(randn(1, V, N), K, 1)); % no diversity

        for pp = 1:V
           for tt = 1:N/2
            x(pp,tt) = w(pp,2*tt) + 2*w(pp,2*tt-1); %ex:10 = 2 01 = 1 %32400
           end
        end
        s = scmaenc(x, CB, h); % joint encoding and fading channel propagation
        y = awgn(s, SNR(k));

        LLR = scmadec(y, CB, h, N0, Niter);
        LLR(LLR==inf) = 1500;
        LLR(LLR==-inf) = -1500;
        % symbol to bit conversion
        %r    = de2bi(x, log2(M), 'left-msb');
        %data = zeros(log2(M)*N/2, V);
        %for kk = 1:V
        %    data(:,kk) = reshape(downsample(r, V, kk-1).',[],1);
        %end

        % LLR to bit conversion
        %datadec = reshape((LLR <= 0), [log2(M) N/2*V]).';
        datar   = zeros(log2(M)*N/2, V);
        for kk = 1:V
            for tt = 1:N/2
                %datar(2*tt-1,kk) = LLR(2*kk-1,tt)*N0*12.5;
                %datar(2*tt,kk) = LLR(2*kk,tt)*N0*12.5;%tanh(rxDataSoft*NOISE_VAR_1D*4.5).^3)*3
                %datar(2*tt-1,kk) = (tanh(LLR(2*kk-1,tt)*N0*2.5).^1)*8;
                %datar(2*tt,kk) = (tanh(LLR(2*kk,tt)*N0*2.5).^1)*8;
                datar(2*tt-1,kk) = LLR(2*kk-1,tt);
                datar(2*tt,kk) = LLR(2*kk,tt);
            end%datar(:,kk) = reshape(downsample(datadec, V, kk-1).', [], 1);
        end
        for pp = 1:V
            %ansbit(pp,:) = nrPolarDecode(datar(:,pp),N/2,N,16,10,false,11);
            ansbit(pp,:) = ldpcDecoder(datar(:,pp));
        end
        err        = sum(xor(dam', ansbit'));
        NBLerr(k) = NBLerr(k) + sum(err>0);
        Nblocks(k) = Nblocks(k) + 6;
        Nerr(:,k)  = Nerr(:,k) + err.';
        Nbits(:,k) = Nbits(:,k) + log2(M)*N*R;
        
    end
    BER(:,k) = Nerr(:,k)./(0.5*Nbits(:,k));
    BLER(k) = NBLerr(k)./Nblocks(k);
    display(sum(BER(:,k))/6)
    
end
plot(EbN0,log10(BER))


function TimeVaryingCIR = gjmi()

    Ln = 1;
    M = 16;
    Ns = 2048;
    Ng = 0;
    CARRIER_FREQ = 3*10^9;
    UE_SPEED = 540;

    sita = zeros(1, Ln);
    JakesAlpha = zeros(Ln, M);
    subAlpha = zeros(1, M);

    for tap = 1:Ln
        sita(tap) = 2 * pi * rand;
        m = 1:M;
        subAlpha = (2 * pi .* m - pi + sita(tap)) ./ (4 * M);
        JakesAlpha(tap, :) = subAlpha;
    end

    initialPhase1 = zeros(1, M);
    initialPhase2 = zeros(1, M);
    for m = 1:M
        initialPhase1(m) = 2 * pi * rand;
        initialPhase2(m) = 2 * pi * rand;
    end

    Fs = Ns*15*1000;
    chip_period = 1/Fs;
    SYMBOL_DURATION = (Ns+Ng)*chip_period;
    dopplerSpread = CARRIER_FREQ * (UE_SPEED * (1000/3600))  / (3* 10^8);
    TimeVaryingCIR = zeros(Ln, (Ns + Ng + Ng));

    % l = 0:Ln-1;
    % tau_function = sqrt((1-0.8)/(1-0.8^Ln))*0.8.^(l./2);
    % for j = 1:length(tau_function)
    %     freq_chann(j) = tau_function(j)*sqrt(0.5)*(randn + 1i * randn);
    % end

    for tap = 1:1
        for idx = 1:(Ns+Ng+Ng)
            timeInstant = (idx-1) * chip_period;
            Inphase = 0; 
            Quadrature = 0;

            for m = 1:M
                Inphase = Inphase + cos(2 * pi * dopplerSpread * timeInstant * cos(JakesAlpha(tap, m)) + initialPhase1(m));
                Quadrature = Quadrature + cos(2 * pi * dopplerSpread * timeInstant * sin(JakesAlpha(tap, m)) + initialPhase2(m));
            end
            TimeVaryingCIR(tap, idx) = (sqrt(1/M) * Inphase + 1i * sqrt(1/M) * Quadrature);
        end
    end
    
    %ppp = ppp + sum(abs(TimeVaryingCIR))/2048

end