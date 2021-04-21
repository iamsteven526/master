clc;
clear;

%ccdf = comm.CCDF('AveragePowerOutputPort',true,'PeakPowerOutputPort',true);
block=10000;

for i = 1: block
    message = randi([0 1],512,1);
    i
    for p = 1:15
        messagenew = message;
        messagenew(1:4) = xor(messagenew(1:4),de2bi(p,4)');
        messagenew = nrPolarEncode(messagenew,1024,10,false);
        qamTxSig = qammod(messagenew,16,'InputType','bit','UnitAveragePower',true);
        TTT = ifft(qamTxSig);
        A(p) = 10*log10(max(abs(TTT).^2)/mean(abs(TTT).^2));
    end
    B(i) = min(A);
end

for j = 1:21
   C(j) = length(B(B>5.75+0.25*j)) ;
end

D = [6:0.25:11];
plot(D,log10(C/block));