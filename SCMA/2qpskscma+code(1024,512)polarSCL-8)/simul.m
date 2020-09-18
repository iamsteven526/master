% Saint Petersburg Electrotechnical University, Saint Petersburg, Russia
% Faculty of Radio Engineering
% Department of Theoretical Fundamentals of Radio Engineering
% Vyacheslav P. Klimentyev and Alexander B. Sergienko, 2015
function BER = simul()
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
          
CB(:,:,7) = [ 0                  0                  0                  0;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071 ];

CB(:,:,8) = [ 0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0 ];

CB(:,:,9) = [0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
              0                  0                  0                  0 ];

CB(:,:,10) = [ 0                  0                  0                  0;...
              0                  0                  0                  0;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071 ];

CB(:,:,11) = [0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0;...
              0                  0                  0                  0;...
             0.7071             1j*0.7071          -1j*0.7071         -0.7071 ];

CB(:,:,12) = [ 0                  0                  0                  0;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0.7071             1j*0.7071          -1j*0.7071         -0.7071;...
              0                  0                  0                  0 ];
         
          

          



K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)

N = 1024; % SCMA signals in frame
R = 0.5;  %code rate
EbN0 = 15:5:25;
SNR  = EbN0 + 10*log10(R*log2(M)*V/K);

Nerr  = zeros(V, length(SNR));
Nbits = zeros(V, length(SNR));
BER   = zeros(V, length(SNR));

maxNumErrs = 10;
maxNumBits = 1e6;
Niter      = 5;
ldpcDecoder = comm.LDPCDecoder;
ldpcEncoder = comm.LDPCEncoder;

for k = 1:length(SNR)

    N0 = 1/(10^(SNR(k)/10)); % noise power

    while ((min(Nerr(:,k)) < maxNumErrs) && (Nbits(1,k) < maxNumBits))
        display(k);
        display(Nbits(1,k));
        %x = randi([0 M-1], V, N); % log2(M)-bit symbols
        dam = randi([0 1], V, N/2); %32400      
        for pp = 1:V
           %w(pp,:) = ldpcEncoder(dam(pp,:)');%64800
           w(pp,:) = nrPolarEncode(dam(pp,:)',N,10,false);
        end
        
        %h = 1/sqrt(2)*(randn(K, V, N)+1j*randn(K, V, N)); % Rayleigh channel
        h = 1/sqrt(2)*(repmat(randn(1, V, N), K, 1)+1j*repmat(randn(1, V, N), K, 1));
        for pp = 1:1
            for qq = 1:N/1
                h(:,:,(pp-1)*N/1+qq) = h(:,:,(pp-1)*N/1+1);
            end
        end
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
        datar = zeros(log2(M)*N/2, V);
        for kk = 1:V
            for tt = 1:N/2
                datar(2*tt-1,kk) = LLR(2*kk-1,tt)*N0/1;
                datar(2*tt,kk) = LLR(2*kk,tt)*N0/1;
            end
        end
        for pp = 1:V
            ansbit(pp,:) = nrPolarDecode(datar(:,pp),N/2,N,8,10,false,24);
            %ansbit(pp,:) = ldpcDecoder(datar(:,pp));
        end

        err        = sum(xor(dam', ansbit'));
        Nerr(:,k)  = Nerr(:,k) + err.';
        Nbits(:,k) = Nbits(:,k) + log2(M)*N*R;
    end
    BER(:,k) = Nerr(:,k)./Nbits(:,k);
    k
    display(sum(Nerr(:,k))/12)
end
plot(EbN0,log10(BER))