clc
clear


bpsk_p1=[1.461712e-01, 6.435860e-02, 2.320917e-02, 7.760750e-03, 2.497190e-03, 8.014200e-04, 2.504000e-04, 8.114000e-05, 2.465000e-05 ];
bpsk_p2=[1.757321e-01, 8.622525e-02, 3.346046e-02, 1.145649e-02, 3.735055e-03, 1.177105e-03, 3.677900e-04, 1.161050e-04, 3.644500e-05 ];
bpsk_p3=[2.039683e-01, 1.150084e-01, 4.851005e-02, 1.730219e-02, 5.717172e-03, 1.834123e-03, 5.786633e-04, 1.837667e-04, 6.015667e-05 ];
bpsk_p4=[2.335615e-01, 1.481179e-01, 7.055732e-02, 2.674499e-02, 9.077484e-03, 2.942479e-03, 9.364925e-04, 2.946125e-04, 9.307250e-05 ];
bpsk_p5=[2.595002e-01, 1.848915e-01, 1.007474e-01, 4.147641e-02, 1.479627e-02, 4.834669e-03, 1.555794e-03, 4.888120e-04, 1.545880e-04 ];
bpsk_p6=[2.823431e-01, 2.212007e-01, 1.379269e-01, 6.395664e-02, 2.397749e-02, 8.064343e-03, 2.593040e-03, 8.225667e-04, 2.608617e-04 ];
bpsk_p7=[3.004057e-01, 2.537012e-01, 1.789330e-01, 9.551838e-02, 3.903127e-02, 1.366024e-02, 4.431711e-03, 1.439708e-03, 4.527943e-04 ];
bpsk_p8=[3.157129e-01, 2.812286e-01, 2.208900e-01, 1.357775e-01, 6.232705e-02, 2.310747e-02, 7.701510e-03, 2.499855e-03, 7.981250e-04 ];
bpsk_p9=[3.273688e-01, 3.034282e-01, 2.576202e-01, 1.805817e-01, 9.553607e-02, 3.882212e-02, 1.353660e-02, 4.408133e-03, 1.410567e-03 ];
bpsk_p10=[3.379101e-01, 3.205122e-01, 2.880826e-01, 2.257239e-01, 1.389225e-01, 6.343419e-02, 2.348428e-02, 7.890900e-03, 2.579920e-03 ];

figure;
semilogy(0:5:40, bpsk_p1, '-o', 'MarkerSize', 10, 'linewidth', 1); hold on;
semilogy(0:5:40, bpsk_p2, '-s', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p3, '-x', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p4, '-+', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p5, '-*', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p6, '-d', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p7, '-<', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p8, '->', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p9, '-p', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_p10, '-h', 'MarkerSize', 10, 'linewidth', 1);

xlabel('$\bar{E_b}/N_0$ [dB]','interpreter','latex', 'fontsize', 20);
ylabel('BER', 'fontsize', 18);
leg = legend('{\itU} = 1; BPSK','{\itU} = 2; BPSK','{\itU} = 3; BPSK','{\itU} = 4; BPSK','{\itU} = 5; BPSK','{\itU} = 6; BPSK','{\itU} = 7; BPSK','{\itU} = 8; BPSK','{\itU} = 9; BPSK','{\itU} = 10; BPSK', 'location', 'southwest');
%leg = legend('{\itP} = 2; SPA','{\itP} = 2; G-SPA','{\itP} = 2; SCL-32','{\itP} = 2; J-SCL-32', 'location', 'southwest');
set(leg, 'fontsize', 15);
set(gca, 'fontname', 'Times New Roman');
set(gcf, 'position', [100, 100, 700, 500]);
grid


bpsk_de_p1=[2.247741e-01, 1.023216e-01, 3.753212e-02, 1.262499e-02, 4.041879e-03, 1.307485e-03, 4.159596e-04, 1.425354e-04, 4.126263e-05];
bpsk_de_p2=[2.616144e-01, 1.372159e-01, 5.533141e-02, 1.929904e-02, 6.339045e-03, 2.017990e-03, 6.359192e-04, 1.967980e-04, 6.961616e-05];
bpsk_de_p3=[2.988254e-01, 1.813619e-01, 8.066861e-02, 3.022297e-02, 9.902140e-03, 3.207057e-03, 1.010902e-03, 3.151448e-04, 1.062593e-04];
bpsk_de_p4=[3.327483e-01, 2.297672e-01, 1.186489e-01, 4.657166e-02, 1.634266e-02, 5.299373e-03, 1.687854e-03, 5.270783e-04, 1.701641e-04];
bpsk_de_p5=[3.615101e-01, 2.797696e-01, 1.665622e-01, 7.379741e-02, 2.666949e-02, 8.837487e-03, 2.836931e-03, 9.112949e-04, 2.880768e-04];
bpsk_de_p6=[3.841208e-01, 3.243098e-01, 2.234595e-01, 1.121599e-01, 4.393932e-02, 1.504573e-02, 4.863601e-03, 1.567135e-03, 4.899040e-04];
bpsk_de_p7=[4.010652e-01, 3.593255e-01, 2.796233e-01, 1.643433e-01, 7.077997e-02, 2.583714e-02, 8.432606e-03, 2.699089e-03, 8.575382e-04];
bpsk_de_p8=[4.147854e-01, 3.863090e-01, 3.297543e-01, 2.257985e-01, 1.120278e-01, 4.363704e-02, 1.483117e-02, 4.802988e-03, 1.538036e-03];

figure;
semilogy(0:5:40, bpsk_de_p1, '-o', 'MarkerSize', 10, 'linewidth', 1); hold on;
semilogy(0:5:40, bpsk_de_p2, '-s', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_de_p3, '-x', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_de_p4, '-+', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_de_p5, '-*', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_de_p6, '-d', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_de_p7, '-<', 'MarkerSize', 10, 'linewidth', 1);
semilogy(0:5:40, bpsk_de_p8, '->', 'MarkerSize', 10, 'linewidth', 1);

xlabel('$\bar{E_b}/N_0$ [dB]','interpreter','latex', 'fontsize', 20);
ylabel('BER', 'fontsize', 18);
leg = legend('{\itU} = 1; DE-BPSK','{\itU} = 2; DE-BPSK','{\itU} = 3; DE-BPSK','{\itU} = 4; DE-BPSK','{\itU} = 5; DE-BPSK','{\itU} = 6; DE-BPSK','{\itU} = 7; DE-BPSK','{\itU} = 8; DE-BPSK', 'location', 'southwest');
%leg = legend('{\itP} = 2; SPA','{\itP} = 2; G-SPA','{\itP} = 2; SCL-32','{\itP} = 2; J-SCL-32', 'location', 'southwest');
set(leg, 'fontsize', 15);
set(gca, 'fontname', 'Times New Roman');
set(gcf, 'position', [100, 100, 700, 500]);
grid