close all;
clear;
clc;

fid=fopen('H_1024x512z32.txt','r');
C = textscan(fid,'%d','delimiter',sprintf(' '));
A = cell2mat(C);
B = reshape(A,1024,512)';
B = sparse(double(B));

M = 2; % Modulation order BPSK)
blk_num = 0;
snr = 3.5;
ldpcEncoder = comm.LDPCEncoder(B);
ldpcDecoder = comm.LDPCDecoder(B);
pskMod = comm.PSKModulator(M,'BitInput',true);
pskDemod = comm.PSKDemodulator(M,'BitOutput',true,...
    'DecisionMethod','Approximate log-likelihood ratio');
pskuDemod = comm.PSKDemodulator(M,'BitOutput',true,...
    'DecisionMethod','Hard decision');
errRate = zeros(1,length(snr));
blkerrRate = zeros(1, length(snr));

for ii = 1:length(snr)
    ttlErr = 0;
    ttlBlkErr = 0;
    pskDemod.Variance = 2/(10^(snr(ii)/10)); % Set variance using current SNR
    while ttlBlkErr < 100 || ttlErr < 1000
        blk_num = blk_num + 1;
        data = logical(randi([0 1],512,1));
        encData = ldpcEncoder(data);
        modSig = pskMod(encData);
        rxSig = awgn(modSig,snr(ii)-3.01,'measured');
        demodSig = pskDemod(rxSig);
        rxBits = ldpcDecoder(demodSig);
        numErr = biterr(data,rxBits);
        if numErr ~= 0
            ttlBlkErr = ttlBlkErr + 1;
        end
        ttlErr = ttlErr + numErr;
    end 

    ttlBits = blk_num*length(rxBits);
    errRate(ii) = ttlErr/ttlBits;
    blkerrRate(ii) = ttlBlkErr/blk_num;
end