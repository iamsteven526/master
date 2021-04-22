clc;
clear;

%ccdf = comm.CCDF('AveragePowerOutputPort',true,'PeakPowerOutputPort',true);
block=100;
Q = 16;
symbolcount = 16;
sccount = 256/symbolcount;

jump = 31;

for i = 1: block
    i
    message = randi([0 1],512,1);
    A = 100;
    for p = 1:2^Q
        messagenew = message;
        messagenew(1:jump:jump*Q) = xor(messagenew(1:jump:jump*Q),de2bi(p-1,Q)');
        messagenew = nrPolarEncode(messagenew,1024,10,false);
        qamTxSig = qammod(messagenew,16,'InputType','bit','UnitAveragePower',true);
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


%origin

%[0.990750000000000,0.977750000000000,0.951250000000000,0.911000000000000,0.862000000000000,0.799750000000000,0.728000000000000,0.640500000000000,0.541750000000000,0.456500000000000,0.381750000000000,0.303500000000000,0.243250000000000,0.178750000000000,0.129250000000000,0.0985000000000000,0.0680000000000000,0.0455000000000000,0.0292500000000000,0.0182500000000000,0.0117500000000000,0.00650000000000000,0.00450000000000000,0.00300000000000000,0.00150000000000000,0.000500000000000000,0.000250000000000000,0,0]
%[0.989087500000000,0.975212500000000,0.951700000000000,0.914300000000000,0.864875000000000,0.799300000000000,0.725300000000000,0.637275000000000,0.549450000000000,0.462762500000000,0.379250000000000,0.301750000000000,0.236387500000000,0.178975000000000,0.133750000000000,0.0968750000000000,0.0684625000000000,0.0470250000000000,0.0310500000000000,0.0198500000000000,0.0123125000000000,0.00715000000000000,0.00426250000000000,0.00228750000000000,0.00118750000000000,0.000562500000000000,0.000237500000000000,6.25000000000000e-05,0]
%Q = 12

%[0.988250000000000,0.970687500000000,0.943437500000000,0.899687500000000,0.840250000000000,0.756125000000000,0.666687500000000,0.555750000000000,0.446437500000000,0.340312500000000,0.246250000000000,0.165687500000000,0.108437500000000,0.0711875000000000,0.0442500000000000,0.0275625000000000,0.0168125000000000,0.00962500000000000,0.00593750000000000,0.00337500000000000,0.00162500000000000,0.000875000000000000,0.000375000000000000,0.000125000000000000,6.25000000000000e-05,6.25000000000000e-05,0,0,0]

%Q=16
%[0.980555555555556,0.961111111111111,0.929861111111111,0.872222222222222,0.815277777777778,0.721527777777778,0.607638888888889,0.462500000000000,0.309722222222222,0.178472222222222,0.106250000000000,0.0576388888888889,0.0354166666666667,0.0215277777777778,0.0125000000000000,0.00555555555555556,0.00347222222222222,0.00138888888888889,0.000694444444444445,0.000694444444444445,0.000694444444444445,0,0,0,0,0,0,0,0]%91
%[0.986111111111111,0.965277777777778,0.930555555555556,0.877976190476191,0.793650793650794,0.689484126984127,0.574404761904762,0.459325396825397,0.320436507936508,0.192460317460317,0.0972222222222222,0.0535714285714286,0.0327380952380952,0.0138888888888889,0.00694444444444444,0.00396825396825397,0.00396825396825397,0.00198412698412698,0.00198412698412698,0.000992063492063492,0,0,0,0,0,0,0,0,0]%63

