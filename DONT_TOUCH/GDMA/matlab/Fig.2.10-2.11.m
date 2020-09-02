clc
clear

pam4_p1=[1.957309e-01, 1.024455e-01, 4.213410e-02, 1.475708e-02, 4.848170e-03, 1.547235e-03, 4.875450e-04, 1.557800e-04, 4.858500e-05];
pam4_p2=[2.400816e-01, 1.539537e-01, 7.857818e-02, 3.223860e-02, 1.146114e-02, 3.770435e-03, 1.216838e-03, 3.747800e-04, 1.145975e-04];
pam4_p3=[2.844804e-01, 2.224614e-01, 1.456594e-01, 7.291407e-02, 2.907387e-02, 1.005263e-02, 3.305043e-03, 1.055637e-03, 3.443000e-04];
pam4_p4=[3.191800e-01, 2.806760e-01, 2.251599e-01, 1.456494e-01, 7.142934e-02, 2.792136e-02, 9.483100e-03, 3.142450e-03, 9.497875e-04];
pam4_p5=[3.448840e-01, 3.220470e-01, 2.898740e-01, 2.325966e-01, 1.492800e-01, 7.170017e-02, 2.764848e-02, 9.483950e-03, 3.053640e-03];

pam8_p1=[2.772831e-01, 1.842674e-01, 1.001075e-01, 4.304792e-02, 1.555465e-02, 5.150050e-03, 1.637123e-03, 5.111433e-04, 1.741067e-04];
pam8_p2=[3.150305e-01, 2.419257e-01, 1.646120e-01, 9.298398e-02, 4.320700e-02, 1.711098e-02, 6.105420e-03, 2.010783e-03, 6.280033e-04];
pam8_p3=[3.499250e-01, 3.085185e-01, 2.689763e-01, 2.150278e-01, 1.397757e-01, 6.906651e-02, 2.742930e-02, 9.514256e-03, 3.092078e-03];


figure;
semilogy(0:5:40, pam4_p1, '-o', 'MarkerSize', 10, 'linewidth', 1); hold on;
semilogy(0:5:40, pam4_p2, '-s', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, pam4_p3, '-x', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, pam4_p4, '-+', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, pam4_p5, '-*', 'MarkerSize', 10, 'linewidth', 1);

xlabel('$\bar{E_b}/N_0$ [dB]','interpreter','latex', 'fontsize', 20);
ylabel('BER', 'fontsize', 18);
leg = legend('{\itU} = 1; 4-PAM','{\itU} = 2; 4-PAM','{\itU} = 3; 4-PAM','{\itU} = 4; 4-PAM','{\itU} = 5; 4-PAM', 'location', 'southwest');
%leg = legend('{\itP} = 2; SPA','{\itP} = 2; G-SPA','{\itP} = 2; SCL-32','{\itP} = 2; J-SCL-32', 'location', 'southwest');
set(leg, 'fontsize', 15);
set(gca, 'fontname', 'Times New Roman');
set(gcf, 'position', [100, 100, 700, 500]);
grid

figure;
semilogy(0:5:40, pam8_p1, '-o', 'MarkerSize', 10, 'linewidth', 1); hold on;
semilogy(0:5:40, pam8_p2, '-s', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, pam8_p3, '-x', 'MarkerSize', 10, 'linewidth', 1);

xlabel('$\bar{E_b}/N_0$ [dB]','interpreter','latex', 'fontsize', 20);
ylabel('BER', 'fontsize', 18);
leg = legend('{\itU} = 1; 8-PAM','{\itU} = 2; 8-PAM','{\itU} = 3; 8-PAM', 'location', 'southwest');
%leg = legend('{\itP} = 2; SPA','{\itP} = 2; G-SPA','{\itP} = 2; SCL-32','{\itP} = 2; J-SCL-32', 'location', 'southwest');
set(leg, 'fontsize', 15);
set(gca, 'fontname', 'Times New Roman');
set(gcf, 'position', [100, 100, 700, 500]);
grid