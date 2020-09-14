clear all;
x=[0:5:40]
user210=[0.204, 0.148, 0.140, 0.134, 0.043,0.0105,0.00354, 0.0011,0.0004];
user220=[0.2343, 0.1368, 0.0608, 0.0231, 0.00920, 0.0023, 0.000702, 0.000328, 0.0000739];
%user210=[0.3022, 0.2343, 0.1368, 0.0756, 0.0125, 0.00281, 0.00103, 0.0003003, 0.0001207];
%user220=[2.822581e-01, 2.205438e-01, 1.383000e-01, 6.415830e-02, 2.411139e-02, 8.022765e-03, 2.555604e-03, 8.369341e-04, 2.718665e-04,];
user230=[1.283102e-01, 4.667958e-02, 1.262283e-02, 4.211750e-03, 1.043000e-03, 3.960417e-04, 1.654167e-04, 4.683333e-05, 0];
user240=[2.044905e-01, 1.149419e-01, 4.821888e-02, 1.734722e-02, 5.721563e-03, 1.755119e-03, 5.961381e-04, 1.831771e-04, 5.632813e-05];
%old = [0.268, 0.184, 0.0987, 0.0377, 0.01286, 0.00454, 0.001233, 0.0004419, 0.0001159];

%bound=0.5*(10.^(-x./10))./(40);
%user210=[1.505403e-01, 3.723727e-02, 8.853785e-03, 2.470017e-03, 9.934857e-04, 5.421531e-04, 4.534110e-04, 4.773200e-04, 4.563461e-04];
%user220=[9.907744e-02, 1.986784e-02, 4.010316e-03, 9.802077e-04, 2.820573e-04, 8.603798e-05, 2.673906e-05, 8.419302e-06, 2.644182e-06];
%user230=[7.735769e-02, 1.416674e-02, 2.681283e-03, 6.503074e-04, 1.836257e-04, 5.588606e-05, 1.734505e-05, 5.522502e-06, 1.735070e-06];
%user240=[6.628320e-02, 1.136525e-02, 2.075260e-03, 4.853708e-04, 1.367947e-04, 4.152985e-05, 1.297536e-05, 4.079959e-06, 1.283563e-06];
figure ;
semilogy(x,user210,'-o', 'MarkerSize', 12, 'linewidth', 1.5);
hold on;
semilogy(x,user220,'-d', 'MarkerSize', 12, 'linewidth', 1.5);
hold on;
semilogy(x,user230,'-<', 'MarkerSize', 12, 'linewidth', 1.5);
hold on;
semilogy(x,user240,'->', 'MarkerSize', 12, 'linewidth', 1.5);
%hold on;
%semilogy(x,old,'->', 'MarkerSize', 12, 'linewidth', 1.5);


xlabel('$\bar{E_b}/N_0$ [dB]', 'interpreter', 'latex', 'fontsize', 23);
ylabel('BER', 'fontsize', 21);

%leg = legend('{\itU} = 3','{\itU} = 4','{\itU} = 5','{\itU} = 6','{\itOld U} = 4','Location','southwest');
leg = legend('QPSK','GDMA U = 3','QPSK with CSI','GDMA U=3 with CSI','Location','southwest')
set(leg, 'fontsize', 15);
set(gca, 'fontname', 'Times New Roman');
set(gcf,'position',[100, 100, 850, 600]);
grid on;
