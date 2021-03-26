clc;
clear;
SNR_dB = 16;% SNR PER BIT
NUM_FRAMES = 10^5; 
FFT_LEN = 512;
NUM_BIT = 4*FFT_LEN; % NUMBER OF DATA BITS
CHAN_LEN = 5; % NUMBER OF CHANNEL TAPS
CP_LEN = CHAN_LEN-1; % LENGTH OF THE CYCLIC PREFIX
FADE_VAR_1D = 1; % 1D FADE VARIANCE OF THE CHANNEL
FADE_STD_DEV = sqrt(FADE_VAR_1D); % STANDARD DEVIATION OF THE FADING CHANNEL
% SNR PER BIT DEFINITION - OVERALL RATE IS 2
SNR = 10^(0.1*SNR_dB); % LINEAR SCALE
R = 2*(256/260);%coderate
NOISE_VAR_1D = 1/(R*2*SNR); % 1D AWGN VARIANCE 
NOISE_STD_DEV = sqrt(NOISE_VAR_1D); % NOISE STANDARD DEVIATION
C_BER = 0; % bit errors in each frame
C_BLER = 0;
st2 = 48511; %random seed
st3 = 87; % random seed

crcLen = 0;
poly = '11';
nPC = 0;
nMax = 10;
iIL = false;
iBIL = true;

interleavepi = 111;


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
    % SOURCE
    A = randi([0 1],1,(NUM_BIT/2) - crcLen);
    %A = nrCRCEncode(A',poly)';
    %enc = nrPolarEncode(A',1024,10,false)';
    enc = ldpcEncoder(A')';
    
    %enc = randintrlv(enc,st2); % Interleave.
    for p = 1:NUM_BIT
       enci(p) = enc(mod(interleavepi*p,NUM_BIT)+1); 
    end
    
    
    % QPSK MAPPING
    F_SIG_NO_CP = qammod(enci',16,'InputType','bit','UnitAveragePower',true)';
    
    %F_SIG_NO_CP = randintrlv(F_SIG_NO_CP,st3);
    % IFFT 
    T_SIG_NO_CP = sqrt(FFT_LEN)*ifft(F_SIG_NO_CP);
    % INSERTING CYCLIC PREFIX
    T_SIG_CP = [T_SIG_NO_CP(end-CP_LEN+1:end) T_SIG_NO_CP];
    %---------------     CHANNEL      -----------------------------------------
    % RAYLEIGH FREQUENCY SELECTIVE FADING CHANNEL
    FADE_CHAN = normrnd(0,FADE_STD_DEV,1,CHAN_LEN)+1i*normrnd(0,FADE_STD_DEV,1,CHAN_LEN);
    normal_power = sqrt(0.5*[0.6364, 0.2341, 0.0861, 0.03168, 0.01165]);
    FADE_CHAN = FADE_CHAN.*normal_power;
    
    FREQ_RESP = fft(FADE_CHAN,FFT_LEN); % ACTUAL CHANNEL FREQUENCY RESPONSE
    % AWGN
    AWGN = normrnd(0,1,1,FFT_LEN+CP_LEN+CHAN_LEN-1)+1i*normrnd(0,1,1,FFT_LEN+CP_LEN+CHAN_LEN-1);
    % CHANNEL OUTPUT
    T_REC_SIG = conv(T_SIG_CP,FADE_CHAN) + NOISE_STD_DEV*AWGN;
    %----------------      RECEIVER  ------------------------------------------
    % CP & TRANSIENT SAMPLES REMOVAL
    T_REC_SIG(1:CP_LEN) = [];
    T_REC_SIG_NO_CP = T_REC_SIG(1:FFT_LEN);
    % PERFORMING THE FFT
    F_REC_SIG_NO_CP = fft(T_REC_SIG_NO_CP)/sqrt(FFT_LEN);
    F_REC_SIG_NO_CP = F_REC_SIG_NO_CP./FREQ_RESP;
    
    %F_REC_SIG_NO_CP = randdeintrlv(F_REC_SIG_NO_CP,st3);
    
    
    %rxDataHard = qamdemod(F_REC_SIG_NO_CP',16,'OutputType','bit','UnitAveragePower',true)';
    rxDataSoft = qamdemod(F_REC_SIG_NO_CP',16,'OutputType','approxllr','UnitAveragePower',true,'NoiseVariance',NOISE_VAR_1D)';
    %rxDataSoft = (rxDataSoft*NOISE_VAR_1D).^3;
    rxDataSoft = (tanh(rxDataSoft*NOISE_VAR_1D*4.5).^3)*3;
    
    for p = 1:NUM_BIT
       rxDataSoftd(mod(interleavepi*p,NUM_BIT)+1) = rxDataSoft(p); 
    end    
    %rxDataSoft = randdeintrlv(rxDataSoft,st2); % Deinterleave.
    
    %rxBits = nrPolarDecode(rxDataSoftd',512,1024,64,10,false,11)';
    rxBits = ldpcDecoder(rxDataSoftd')';
    
    numErrsInFrameHard = biterr(A(1:(NUM_BIT/2)-crcLen),rxBits(1:(NUM_BIT/2)-crcLen));
    if numErrsInFrameHard >= 1
        bler_flag = 1;
    else
        bler_flag = 0;
    end
    C_BER = C_BER + numErrsInFrameHard;
    C_BLER = C_BLER + bler_flag;
end
BER = C_BER/(NUM_BIT*NUM_FRAMES);
BLER = C_BLER/NUM_FRAMES;


