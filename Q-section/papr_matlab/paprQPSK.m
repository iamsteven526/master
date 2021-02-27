clc;
clear;

%ccdf = comm.CCDF('AveragePowerOutputPort',true,'PeakPowerOutputPort',true);
block=100;
Q = 16;
symbolcount = 8;
sccount = 512/symbolcount;

jump = 31;

for i = 1: block
    message = randi([0 1],512,1);
    A = 100;
    for p = 1:2^Q
        messagenew = message;
        messagenew(1:jump:jump*Q) = xor(messagenew(1:jump:jump*Q),de2bi(p-1,Q)');
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

%32sym origin:
%[0.991900000000000,0.981100000000000,0.951200000000000,0.916100000000000,0.848100000000000,0.795400000000000,0.739200000000000,0.669100000000000,0.554300000000000,0.445700000000000,0.391100000000000,0.312600000000000,0.261500000000000,0.209300000000000,0.151500000000000,0.112300000000000,0.0801000000000000,0.0559000000000000,0.0391000000000000,0.0269000000000000,0.0164000000000000,0.00990000000000000,0.00740000000000000,0.00430000000000000,0.00280000000000000,0.00130000000000000,0.000400000000000000,0.000200000000000000,0.000200000000000000]





%Q=12,32sym:
%[0.991593750000000,0.982343750000000,0.950656250000000,0.913812500000000,0.843031250000000,0.788406250000000,0.733937500000000,0.668750000000000,0.555156250000000,0.445937500000000,0.390812500000000,0.308375000000000,0.253500000000000,0.201468750000000,0.144625000000000,0.100531250000000,0.0696562500000000,0.0473750000000000,0.0337500000000000,0.0217812500000000,0.0132500000000000,0.00731250000000000,0.00471875000000000,0.00287500000000000,0.00196875000000000,0.00112500000000000,0.000625000000000000,0.000218750000000000,9.37500000000000e-05]

%Q=12,32sym:
%[0.992187500000000,0.981875000000000,0.950156250000000,0.913203125000000,0.838593750000000,0.778984375000000,0.723437500000000,0.654921875000000,0.524609375000000,0.407656250000000,0.344140625000000,0.250156250000000,0.187265625000000,0.131015625000000,0.0775781250000000,0.0437500000000000,0.0244531250000000,0.0139843750000000,0.00921875000000000,0.00484375000000000,0.00234375000000000,0.00101562500000000,0.000625000000000000,7.81250000000000e-05,0,0,0,0,0]

%Q=16,32sym:
%[0.992127862595420,0.982347328244275,0.944895038167939,0.906011450381679,0.824666030534351,0.764551526717557,0.701812977099237,0.631679389312977,0.492604961832061,0.369751908396947,0.308206106870229,0.207776717557252,0.142891221374046,0.0839694656488550,0.0391221374045802,0.0214694656488550,0.0104961832061069,0.00572519083969466,0.00357824427480916,0.00214694656488550,0.000477099236641221,0.000238549618320611,0.000238549618320611,0,0,0,0,0,0]





%Q=8,32sym:
%[0.991687500000000,0.980187500000000,0.948312500000000,0.911062500000000,0.839593750000000,0.786125000000000,0.734375000000000,0.668468750000000,0.555093750000000,0.447781250000000,0.393437500000000,0.306093750000000,0.252937500000000,0.200625000000000,0.143187500000000,0.102593750000000,0.0705312500000000,0.0488125000000000,0.0342187500000000,0.0221250000000000,0.0129375000000000,0.00709375000000000,0.00493750000000000,0.00275000000000000,0.00184375000000000,0.00106250000000000,0.000531250000000000,0.000218750000000000,9.37500000000000e-05]


%16sym origin
%16sym
%[0.999950000000000,0.999850000000000,0.998600000000000,0.994900000000000,0.986200000000000,0.970350000000000,0.938550000000000,0.898600000000000,0.833650000000000,0.745700000000000,0.671800000000000,0.583150000000000,0.494700000000000,0.397850000000000,0.317850000000000,0.242850000000000,0.186350000000000,0.137200000000000,0.0960500000000000,0.0691000000000000,0.0455500000000000,0.0293000000000000,0.0187500000000000,0.0117000000000000,0.00660000000000000,0.00400000000000000,0.00250000000000000,0.00145000000000000,0.000600000000000000]

%Q=8,16sym:
%[0.999937500000000,0.999687500000000,0.998187500000000,0.994062500000000,0.984812500000000,0.968062500000000,0.937187500000000,0.898625000000000,0.831937500000000,0.743875000000000,0.664750000000000,0.571125000000000,0.479000000000000,0.379312500000000,0.291312500000000,0.216875000000000,0.160250000000000,0.117937500000000,0.0791875000000000,0.0552500000000000,0.0356875000000000,0.0222500000000000,0.0141875000000000,0.00912500000000000,0.00575000000000000,0.00362500000000000,0.00256250000000000,0.00156250000000000,0.00100000000000000]

%Q=12,16sym:
%[0.999921679197995,0.999686716791980,0.998355263157895,0.993969298245614,0.982534461152882,0.962640977443609,0.926456766917293,0.877662907268170,0.795817669172932,0.689771303258145,0.599232456140351,0.483552631578947,0.368264411027569,0.244595864661654,0.155231829573935,0.0890507518796992,0.0527882205513785,0.0303884711779449,0.0180921052631579,0.0111998746867168,0.00681390977443609,0.00407268170426065,0.00258458646616541,0.00133145363408521,0.000783208020050125,0.000548245614035088,0.000313283208020050,0.000234962406015038,7.83208020050125e-05]



%origin8sym
%[1,1,1,1,0.999900000000000,0.999200000000000,0.997450000000000,0.991000000000000,0.974200000000000,0.943950000000000,0.898700000000000,0.838250000000000,0.753550000000000,0.653900000000000,0.550900000000000,0.446700000000000,0.356100000000000,0.274250000000000,0.201750000000000,0.148850000000000,0.103700000000000,0.0694500000000000,0.0466000000000000,0.0322000000000000,0.0190000000000000,0.0107500000000000,0.00540000000000000,0.00270000000000000,0.00110000000000000]
%Q=8,8sym:
%[1,1,1,1,0.999875000000000,0.999500000000000,0.997000000000000,0.991125000000000,0.974750000000000,0.942250000000000,0.895750000000000,0.828000000000000,0.744250000000000,0.626250000000000,0.520000000000000,0.410125000000000,0.319375000000000,0.237875000000000,0.168500000000000,0.116250000000000,0.0781250000000000,0.0543750000000000,0.0353750000000000,0.0216250000000000,0.0143750000000000,0.00875000000000000,0.00487500000000000,0.00250000000000000,0.00162500000000000]

%[1,1,1,1,0.999875000000000,0.999250000000000,0.996625000000000,0.991000000000000,0.974375000000000,0.941000000000000,0.894750000000000,0.828750000000000,0.738687500000000,0.627062500000000,0.516875000000000,0.405250000000000,0.307875000000000,0.222562500000000,0.152187500000000,0.106562500000000,0.0709375000000000,0.0460000000000000,0.0301250000000000,0.0195625000000000,0.0115625000000000,0.00775000000000000,0.00493750000000000,0.00318750000000000,0.00156250000000000]

%jump=37,Q=8:
%[1,1,1,1,0.999928977272727,0.998721590909091,0.996235795454545,0.988707386363636,0.972869318181818,0.934943181818182,0.877059659090909,0.795454545454545,0.690411931818182,0.560511363636364,0.429545454545455,0.299928977272727,0.204261363636364,0.137286931818182,0.0857954545454546,0.0570312500000000,0.0355823863636364,0.0212357954545455,0.0129261363636364,0.00774147727272727,0.00490056818181818,0.00262784090909091,0.00113636363636364,0.000497159090909091,0.000142045454545455]

%jump=37,Q=12:
%[1,1,1,0.999968750000000,0.999750000000000,0.998531250000000,0.995250000000000,0.985750000000000,0.961093750000000,0.907968750000000,0.829593750000000,0.712250000000000,0.561593750000000,0.373656250000000,0.217500000000000,0.114531250000000,0.0596250000000000,0.0309375000000000,0.0156875000000000,0.00825000000000000,0.00450000000000000,0.00218750000000000,0.00103125000000000,0.000562500000000000,0.000218750000000000,0.000125000000000000,3.12500000000000e-05,3.12500000000000e-05,0]

%jump=31,Q=16:
%[1,1,1,1,1,0.998062015503876,0.992248062015504,0.981589147286822,0.947674418604651,0.884689922480620,0.775193798449612,0.622093023255814,0.391472868217054,0.172480620155039,0.0843023255813954,0.0436046511627907,0.0232558139534884,0.0155038759689922,0.0116279069767442,0.00968992248062016,0.00581395348837209,0.00387596899224806,0.00387596899224806,0.00193798449612403,0,0,0,0,0]
%[1,1,1,1,1,0.997950819672131,0.995901639344262,0.973360655737705,0.940573770491803,0.870901639344262,0.756147540983607,0.614754098360656,0.389344262295082,0.176229508196721,0.0635245901639344,0.0225409836065574,0.0122950819672131,0.00819672131147541,0.00204918032786885,0,0,0,0,0,0,0,0,0,0]



%Q=8,4sym:
%[1,1,1,1,1,1,1,0.999875000000000,0.999625000000000,0.997875000000000,0.990000000000000,0.973625000000000,0.938875000000000,0.874750000000000,0.784625000000000,0.682125000000000,0.557875000000000,0.438000000000000,0.325000000000000,0.236625000000000,0.167375000000000,0.114500000000000,0.0767500000000000,0.0532500000000000,0.0347500000000000,0.0210000000000000,0.0138750000000000,0.00612500000000000,0.00362500000000000]

%Q=8,2sym:
%[1,1,1,1,1,1,1,1,1,1,1,0.998250000000000,0.989500000000000,0.959000000000000,0.901750000000000,0.810250000000000,0.695750000000000,0.558750000000000,0.433000000000000,0.306750000000000,0.224250000000000,0.157750000000000,0.106500000000000,0.0672500000000000,0.0415000000000000,0.0265000000000000,0.0167500000000000,0.00875000000000000,0.00475000000000000]


%Q=4,2sym:
%[1,1,1,1,1,1,1,1,1,1,0.999750000000000,0.998500000000000,0.994000000000000,0.977500000000000,0.929250000000000,0.833250000000000,0.705500000000000,0.566500000000000,0.436750000000000,0.325500000000000,0.227500000000000,0.161750000000000,0.105000000000000,0.0665000000000000,0.0412500000000000,0.0250000000000000,0.0165000000000000,0.0100000000000000,0.00600000000000000]

%Q=4,1sym:
%[1,1,1,1,1,1,1,1,1,1,1,1,0.999500000000000,0.999000000000000,0.978000000000000,0.906000000000000,0.721500000000000,0.469000000000000,0.227000000000000,0.100500000000000,0.0390000000000000,0.0185000000000000,0.00700000000000000,0.00350000000000000,0.00100000000000000,0,0,0,0]



