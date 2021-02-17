clc;
clear;

%ccdf = comm.CCDF('AveragePowerOutputPort',true,'PeakPowerOutputPort',true);
block=1000;
Q = 2;
symbolcount = 1;
sccount = 512/symbolcount;

for i = 1: block
    message = randi([0 1],512,1);
    A = 100;
    for p = 1:2^Q
        messagenew = message;
        messagenew(1:Q) = xor(messagenew(1:Q),de2bi(p-1,Q)');
        messagenew = nrPolarEncode(messagenew,1024,10,false);
        qamTxSig = qammod(messagenew,4,'InputType','bit','UnitAveragePower',true);
        for sym = 1:symbolcount
            TTT = ifft(qamTxSig(sccount*(sym-1)+1:sccount*(sym-1)+sccount));
            G(sym) = 10*log10(max(abs(TTT).^2)/mean(abs(TTT).^2));
        end
        
        if max(G) < max(A)
           A = G; 
        end
        
    end
    B((i-1)*symbolcount+1:(i-1)*symbolcount+symbolcount) = A;
end

for j = 1:29
   C(j) = length(B(B>2.75+0.25*j)) ;
end

D = [3:0.25:10];
test = C/(block*symbolcount);
plot(D,log10(C/(block*symbolcount)));