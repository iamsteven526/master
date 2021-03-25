clc; clear;
Ln = 1;
M = 16;
Ns = 1024;
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

channelMatrix = zeros(Ns, Ns);
for time = 0:Ns-1
    for l = 0:Ln-1
        channelMatrix(time+1, mod(time - l + Ns, Ns) + 1) = TimeVaryingCIR(l+1, time + 1);
    end
end
%H = fft(freq_chann, Ns);
H = (1/Ns) * dftmtx(Ns) * channelMatrix * dftmtx(Ns)';

figure; hold on;
%plot(real(H), imag(H), 'LineWidth',2);
x = (1:Ns+Ng+Ng);
plot(x,abs(TimeVaryingCIR), 'LineWidth',2);
grid on;    box on ; hold off;