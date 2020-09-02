clc
clear

qpsk_p1=[1.461795e-01, 6.422421e-02, 2.322346e-02, 7.749505e-03, 2.527930e-03, 7.704550e-04, 2.465000e-04, 8.095500e-05, 2.308500e-05];
qpsk_p2=[2.203803e-01, 1.321858e-01, 5.951980e-02, 2.206015e-02, 7.481600e-03, 2.294725e-03, 7.458750e-04, 2.220000e-04, 7.162500e-05];
qpsk_p3=[2.767657e-01, 2.133108e-01, 1.290956e-01, 5.809593e-02, 2.151675e-02, 7.072550e-03, 2.342550e-03, 7.262667e-04, 2.328833e-04];
qpsk_p4=[3.136304e-01, 2.780857e-01, 2.159207e-01, 1.304638e-01, 5.919003e-02, 2.179136e-02, 7.327775e-03, 2.342713e-03, 7.358750e-04];
qpsk_p5=[3.375435e-01, 3.197179e-01, 2.862355e-01, 2.222234e-01, 1.351267e-01, 6.151032e-02, 2.290262e-02, 7.660700e-03, 2.479790e-03];

figure;
semilogy(0:5:40, qpsk_p1, '-o', 'MarkerSize', 10, 'linewidth', 1); hold on;
semilogy(0:5:40, qpsk_p2, '-s', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, qpsk_p3, '-x', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, qpsk_p4, '-+', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, qpsk_p5, '-*', 'MarkerSize', 10, 'linewidth', 1);

xlabel('$\bar{E_b}/N_0$ [dB]','interpreter','latex', 'fontsize', 20);
ylabel('BER', 'fontsize', 18);
leg = legend('{\itU} = 1; QPSK','{\itU} = 2; QPSK','{\itU} = 3; QPSK','{\itU} = 4; QPSK','{\itU} = 5; QPSK', 'location', 'southwest');
%leg = legend('{\itP} = 2; SPA','{\itP} = 2; G-SPA','{\itP} = 2; SCL-32','{\itP} = 2; J-SCL-32', 'location', 'southwest');
set(leg, 'fontsize', 15);
set(gca, 'fontname', 'Times New Roman');
set(gcf, 'position', [100, 100, 700, 500]);
grid