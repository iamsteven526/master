clc;
clear;
SNR_dB = 14;% SNR PER BIT (EbN0)
NUM_FRAMES = 12050; 

FFT_LEN = 512*4;
NUM_BIT = 512; % NUMBER OF DATA BITS
CHAN_LEN = 5; % NUMBER OF CHANNEL TAPS
CP_LEN = CHAN_LEN-1; % LENGTH OF THE CYCLIC PREFIX
FADE_VAR_1D = 1; % 1D FADE VARIANCE OF THE CHANNEL
FADE_STD_DEV = sqrt(FADE_VAR_1D); % STANDARD DEVIATION OF THE FADING CHANNEL
% SNR PER BIT DEFINITION - OVERALL RATE IS 2


SNR = 10^(0.1*SNR_dB); % LINEAR SCALE
R = (16/20)*(501/512);%coderate
NOISE_VAR_1D = 1/(R*SNR); % 1D AWGN VARIANCE 
NOISE_STD_DEV = sqrt(NOISE_VAR_1D); % NOISE STANDARD DEVIATION
C_BER = 0; % bit errors in each frame
C_BLER = 0;
st2 = 48511; %random seed
st3 = 87; % random seed

crcLen = 11;
poly = '11';
nPC = 0;
nMax = 10;
iIL = false;
iBIL = true;

interleavepi = 111;


%codebook
CB(:,:,1) = [ 0                  0                  0                  0;...
             (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
              0                  0                  0                  0;...
              (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j)];

CB(:,:,2) = [(-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
              0                  0                  0                  0;...
             (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
              0                  0                  0                  0 ];

CB(:,:,3) = [(-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
             (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
              0                  0                  0                  0;...
              0                  0                  0                  0 ];

CB(:,:,4) = [ 0                  0                  0                  0;...
              0                  0                  0                  0;...
              (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
             (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j) ];

CB(:,:,5) = [(-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
              0                  0                  0                  0;...
              0                  0                  0                  0;...
             (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j)];

CB(:,:,6) = [ 0                  0                  0                  0;...
              (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
              (-0.7071 - 0.7071*1j) (-0.7071 + 0.7071*1j) (0.7071 - 0.7071*1j)      (0.7071 + 0.7071*1j);...
              0                  0                  0                  0 ];

%CB=CB*0.8163*0.7071;         


K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)
%endofcodebook



fid=fopen('H_2048_1024_z64_0635.txt','r');
C = textscan(fid,'%d','delimiter',sprintf(' '));
A = cell2mat(C);
B = reshape(A,2048,1024)';
B = sparse(double(B));

ldpcDecoder = comm.LDPCDecoder(B);
ldpcDecoder.MaximumIterationCount = 250;
ldpcDecoder.IterationTerminationCondition = 'Parity check satisfied';
ldpcEncoder = comm.LDPCEncoder(B);


for FRAME_CNT = 1:NUM_FRAMES
    FRAME_CNT
    % SOURCE
    A = randi([0 1],V,(NUM_BIT) - crcLen);
    A = nrCRCEncode(A',poly)';

    for pp = 1:V
       enc(pp,:) = nrPolarEncode(A(pp,:)',1024,10,false)';
       enc(pp,:) = randintrlv(enc(pp,:),st2);
       %enc(pp,:) = ldpcEncoder(A(pp,:)');
    end    
    
    % QPSK MAPPING
    F_SIG_NO_CP_OLD = qammod(enc',4,'InputType','bit','UnitAveragePower',true)';
    
    F_SIG_NO_CP(1,:) = [zeros(1,NUM_BIT) F_SIG_NO_CP_OLD(1,:) zeros(1,NUM_BIT) F_SIG_NO_CP_OLD(1,:)];
    F_SIG_NO_CP(2,:) = [F_SIG_NO_CP_OLD(2,:) zeros(1,NUM_BIT) F_SIG_NO_CP_OLD(2,:) zeros(1,NUM_BIT)];
    F_SIG_NO_CP(3,:) = [F_SIG_NO_CP_OLD(3,:) F_SIG_NO_CP_OLD(3,:) zeros(1,NUM_BIT) zeros(1,NUM_BIT)];
    F_SIG_NO_CP(4,:) = [zeros(1,NUM_BIT) zeros(1,NUM_BIT) F_SIG_NO_CP_OLD(4,:) F_SIG_NO_CP_OLD(4,:)];
    F_SIG_NO_CP(5,:) = [F_SIG_NO_CP_OLD(5,:) zeros(1,NUM_BIT) zeros(1,NUM_BIT) F_SIG_NO_CP_OLD(5,:)];
    F_SIG_NO_CP(6,:) = [zeros(1,NUM_BIT) F_SIG_NO_CP_OLD(6,:) F_SIG_NO_CP_OLD(6,:) zeros(1,NUM_BIT)];
    
    % IFFT
    T_REC_SIG = zeros(128,16+CP_LEN+CHAN_LEN-1);
    for pp = 1:V
        FADE_CHAN(pp,:) = normrnd(0,FADE_STD_DEV,1,CHAN_LEN)+1i*normrnd(0,FADE_STD_DEV,1,CHAN_LEN);
        normal_power = sqrt(0.5*[0.6364, 0.2341, 0.0861, 0.03168, 0.01165]);
        FADE_CHAN(pp,:) = FADE_CHAN(pp,:).*normal_power;
        FREQ_RESP(pp,:) = fft(FADE_CHAN(pp,:),16);
        for symnum = 1:128
            T_SIG_NO_CP(pp,:) = sqrt(16)*ifft(F_SIG_NO_CP(pp,symnum:128:2048));
            T_SIG_CP(pp,:) = [T_SIG_NO_CP(pp,end-CP_LEN+1:end) T_SIG_NO_CP(pp,:)];

            T_REC_SIG(symnum,:) = T_REC_SIG(symnum,:) + conv(T_SIG_CP(pp,:),FADE_CHAN(pp,:));
        end      
    end
    % INSERTING CYCLIC PREFIX
    
    %---------------     CHANNEL      -----------------------------------------
    % RAYLEIGH FREQUENCY SELECTIVE FADING CHANNEL
    %FADE_CHAN = normrnd(0,FADE_STD_DEV,1,CHAN_LEN)+1i*normrnd(0,FADE_STD_DEV,1,CHAN_LEN);
    %normal_power = sqrt(0.5*[0.6364, 0.2341, 0.0861, 0.03168, 0.01165]);
    %FADE_CHAN = FADE_CHAN.*normal_power;
    
    %FREQ_RESP = fft(FADE_CHAN,FFT_LEN); % ACTUAL CHANNEL FREQUENCY RESPONSE
    % AWGN
    AWGN = normrnd(0,1,128,16+CP_LEN+CHAN_LEN-1)+1i*normrnd(0,1,128,16+CP_LEN+CHAN_LEN-1);
    % CHANNEL OUTPUT
    T_REC_SIG =  T_REC_SIG + NOISE_STD_DEV*AWGN;
    
    
    
    
    %----------------      RECEIVER  ------------------------------------------
    % CP & TRANSIENT SAMPLES REMOVAL
    T_REC_SIG(:,1:CP_LEN) = [];
    for symnum = 1:128
        T_REC_SIG_NO_CP(symnum,:) = T_REC_SIG(symnum,1:16);
        % PERFORMING THE FFT
        F_REC_SIG_NO_CP(symnum,:) = fft(T_REC_SIG_NO_CP(symnum,:))/sqrt(16);
    end
    %TODO: modify input of scmadec() y , CB , h
    h = zeros(K, V, NUM_BIT); % Rayleigh channel
    F_REC_SIG_NO_CP_P = reshape(F_REC_SIG_NO_CP,1,2048);
    for kkk = 1:K
        for vvv = 1:V
            for sc = 1:4
                h(kkk,vvv,(sc-1)*128+1:128*sc) = FREQ_RESP(vvv,4*(kkk-1)+sc);
            end
        end
    end
    y = zeros(K,NUM_BIT);
    y(1,:) = F_REC_SIG_NO_CP_P(1:NUM_BIT);
    y(2,:) = F_REC_SIG_NO_CP_P(NUM_BIT+1:NUM_BIT*2);
    y(3,:) = F_REC_SIG_NO_CP_P(2*NUM_BIT+1:NUM_BIT*3);
    y(4,:) = F_REC_SIG_NO_CP_P(3*NUM_BIT+1:NUM_BIT*4);
    
    LLR = scmadec(y, CB, h, NOISE_STD_DEV*NOISE_STD_DEV , 5);
    LLR(LLR==inf) = 1500;
    LLR(LLR==-inf) = -1500;   
    
    datar   = zeros(2*NUM_BIT, V);
    for kk = 1:V
        for tt = 1:NUM_BIT
            datar(2*tt-1,kk) = LLR(2*kk-1,tt);%(tanh(LLR(2*kk-1,tt)*N0*2.5).^1)*8;
            datar(2*tt,kk) = LLR(2*kk,tt);%(tanh(LLR(2*kk,tt)*N0*2.5).^1)*8;
        end
    end    
    
    
    
    for pp = 1:V
        %ansbit(pp,:) = ldpcDecoder(datar(:,pp));
        datar(:,pp) = randdeintrlv(datar(:,pp),st2);
        ansbit(pp,:) = nrPolarDecode(datar(:,pp),512,1024,32,10,false,11);
    end  
    
    err        = sum(xor(A(:,1:512-crcLen)', double(ansbit(:,1:512-crcLen)')));    
    
%     if err >= 1
%         bler_flag = 1;
%     else
%         bler_flag = 0;
%     end
    bler_flag = sum(err>0);
    C_BER = C_BER + sum(err);
     C_BLER = C_BLER + bler_flag;
end
BER = C_BER/(NUM_BIT*6*NUM_FRAMES);
BLER = C_BLER/(6*NUM_FRAMES);
