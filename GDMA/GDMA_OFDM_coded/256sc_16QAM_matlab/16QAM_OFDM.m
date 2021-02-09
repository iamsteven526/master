clc;
clear;
SNR_dB = 40;% SNR PER BIT
NUM_FRAMES = 10^4; 
FFT_LEN = 256;
NUM_BIT = 4*FFT_LEN; % NUMBER OF DATA BITS
CHAN_LEN = 5; % NUMBER OF CHANNEL TAPS
CP_LEN = CHAN_LEN-1; % LENGTH OF THE CYCLIC PREFIX
FADE_VAR_1D = 0.5; % 1D FADE VARIANCE OF THE CHANNEL
FADE_STD_DEV = sqrt(FADE_VAR_1D); % STANDARD DEVIATION OF THE FADING CHANNEL
% SNR PER BIT DEFINITION - OVERALL RATE IS 2
SNR = 10^(0.1*SNR_dB); % LINEAR SCALE
NOISE_VAR_1D = 0.5*2*2*CHAN_LEN*FADE_VAR_1D/(2*SNR*FFT_LEN); % 1D AWGN VARIANCE 
NOISE_STD_DEV = sqrt(NOISE_VAR_1D); % NOISE STANDARD DEVIATION
C_BER = 0; % bit errors in each frame
for FRAME_CNT = 1:NUM_FRAMES
    % SOURCE
    A = randi([0 1],1,NUM_BIT/2);
    enc = nrPolarEncode(A',1024,10,false)';
    % QPSK MAPPING
    F_SIG_NO_CP = qammod(enc',16,'InputType','bit','UnitAveragePower',true)';
    % IFFT 
    T_SIG_NO_CP = ifft(F_SIG_NO_CP);
    % INSERTING CYCLIC PREFIX
    T_SIG_CP = [T_SIG_NO_CP(end-CP_LEN+1:end) T_SIG_NO_CP];
    %---------------     CHANNEL      -----------------------------------------
    % RAYLEIGH FREQUENCY SELECTIVE FADING CHANNEL
    FADE_CHAN = normrnd(0,FADE_STD_DEV,1,CHAN_LEN)+1i*normrnd(0,FADE_STD_DEV,1,CHAN_LEN);
    normal_power = [0.6364, 0.2341, 0.0861, 0.03168, 0.01165];
    FADE_CHAN = FADE_CHAN.*normal_power;
    
    FREQ_RESP = fft(FADE_CHAN,FFT_LEN); % ACTUAL CHANNEL FREQUENCY RESPONSE
    % AWGN
    AWGN = normrnd(0,NOISE_STD_DEV,1,FFT_LEN+CP_LEN+CHAN_LEN-1)+1i*normrnd(0,NOISE_STD_DEV,1,FFT_LEN+CP_LEN+CHAN_LEN-1);
    % CHANNEL OUTPUT
    T_REC_SIG = conv(T_SIG_CP,FADE_CHAN) + AWGN;
    %----------------      RECEIVER  ------------------------------------------
    % CP & TRANSIENT SAMPLES REMOVAL
    T_REC_SIG(1:CP_LEN) = [];
    T_REC_SIG_NO_CP = T_REC_SIG(1:FFT_LEN);
    % PERFORMING THE FFT
    F_REC_SIG_NO_CP = fft(T_REC_SIG_NO_CP);
    F_REC_SIG_NO_CP = F_REC_SIG_NO_CP./FREQ_RESP;

    %rxDataHard = qamdemod(F_REC_SIG_NO_CP',16,'OutputType','bit','UnitAveragePower',true)';
    rxDataSoft = qamdemod(F_REC_SIG_NO_CP',16,'OutputType','approxllr','UnitAveragePower',true,'NoiseVariance',NOISE_VAR_1D);
    rxBits = nrPolarDecode(rxDataSoft,512,1024,32,10,false,11)';
    
    
    numErrsInFrameHard = biterr(A,rxBits);
    C_BER = C_BER + numErrsInFrameHard;
end
BER = C_BER/(NUM_BIT*NUM_FRAMES);


