clc
clear

G=[0:0.1:3]
S_slotted_ALOHA=[0,9.059100e-02, 1.635130e-01, 2.223770e-01, 2.684039e-01, 3.041636e-01, 3.298512e-01, 3.468318e-01, 3.594371e-01, 3.637922e-01, 3.684519e-01, 3.640812e-01, 3.597225e-01, 3.553598e-01, 3.453515e-01, 3.327907e-01, 3.236421e-01, 3.110139e-01, 2.967316e-01, 2.861327e-01, 2.716897e-01, 2.571507e-01, 2.430664e-01, 2.312863e-01, 2.173799e-01, 2.043836e-01, 1.918990e-01, 1.821608e-01, 1.716230e-01, 1.593151e-01, 1.483429e-01];
S_slotted_ALOHA_GDMA=[0,9.710000e-02, 2.023000e-01, 3.010000e-01, 3.950000e-01, 5.079000e-01, 6.016000e-01, 7.141000e-01, 8.019000e-01, 8.918000e-01, 9.863000e-01, 1.101600e+00, 1.194900e+00, 1.290800e+00, 1.398600e+00, 1.508900e+00, 1.610300e+00, 1.678900e+00, 1.764800e+00, 1.902000e+00, 1.971100e+00, 2.098600e+00, 2.165100e+00, 2.305100e+00, 2.363900e+00, 2.489900e+00, 2.547600e+00, 2.654200e+00, 2.751500e+00, 2.806700e+00, 2.901000e+00];
S_slotted_ALOHA_GDMA_es=[0,1.045000e-01, 1.983000e-01, 2.996000e-01, 4.030000e-01, 4.953000e-01, 5.994000e-01, 6.784000e-01, 7.657000e-01, 8.700000e-01, 9.300000e-01, 9.825000e-01, 1.066900e+00, 1.147900e+00, 1.193100e+00, 1.253000e+00, 1.299400e+00, 1.320600e+00, 1.379700e+00, 1.401100e+00, 1.425300e+00, 1.446400e+00, 1.471500e+00, 1.503200e+00, 1.516500e+00, 1.503200e+00, 1.522800e+00, 1.519800e+00, 1.501600e+00, 1.516600e+00, 1.503100e+00];
S_unslotted_ALOHA=[0,8.086500e-02, 1.329850e-01, 1.638611e-01, 1.781842e-01, 1.826948e-01, 1.801816e-01, 1.719636e-01, 1.613680e-01, 1.503855e-01, 1.345013e-01, 1.224371e-01, 1.096552e-01, 9.668667e-02, 8.542484e-02, 7.582511e-02, 6.516768e-02, 5.995130e-02, 4.863874e-02, 4.361468e-02, 3.803119e-02, 3.326531e-02, 2.769892e-02, 2.483221e-02, 2.035545e-02, 1.810811e-02, 1.594388e-02, 1.335979e-02, 1.126374e-02, 9.197708e-03, 8.727811e-03];
%S_unslotted_ALOHA_GDMA=[0,8.474000e-02, 1.549900e-01, 2.154300e-01, 2.679800e-01, 3.112700e-01, 3.441900e-01, 3.734400e-01, 3.965800e-01, 4.134300e-01, 4.259900e-01, 4.331500e-01, 4.373000e-01, 4.405800e-01, 4.378700e-01, 4.365200e-01, 4.282700e-01, 4.228400e-01, 4.102800e-01, 4.013100e-01, 3.932800e-01, 3.789300e-01, 3.680500e-01, 3.595300e-01, 3.478100e-01, 3.321200e-01, 3.210700e-01, 3.080700e-01, 2.959200e-01, 2.811000e-01, 2.692000e-01];
S_unslotted_ALOHA_GDMA=[0,8.187181e-02, 1.570443e-01, 2.229877e-01, 2.831017e-01, 3.340166e-01, 3.832617e-01, 4.235476e-01, 4.588741e-01, 4.879412e-01, 5.153785e-01, 5.402060e-01, 5.576842e-01, 5.714129e-01, 5.873913e-01, 5.929807e-01, 6.008799e-01, 6.041896e-01, 6.088491e-01, 6.055194e-01, 6.045195e-01, 6.063994e-01, 6.013199e-01, 5.919708e-01, 5.852215e-01, 5.822518e-01, 5.714129e-01, 5.655134e-01, 5.562744e-01, 5.449755e-01, 5.349965e-01];
%S_unslotted_time_estimation=[0,8.170000e-02, 1.531000e-01, 2.142700e-01, 2.631200e-01, 3.005800e-01, 3.376100e-01, 3.623600e-01, 3.851400e-01, 4.002400e-01, 4.098200e-01, 4.158700e-01, 4.176800e-01, 4.198100e-01, 4.161000e-01, 4.116200e-01, 4.076000e-01, 3.986900e-01, 3.906500e-01, 3.808200e-01, 3.662700e-01, 3.600700e-01, 3.436200e-01, 3.311500e-01, 3.193900e-01, 3.119800e-01, 3.002700e-01, 2.851100e-01, 2.733200e-01, 2.623800e-01, 2.528500e-01];
%S_unslotted_time_estimation=[0,7.926207e-02, 1.470253e-01, 2.024398e-01, 2.484352e-01, 2.813119e-01, 3.142986e-01, 3.359364e-01, 3.523648e-01, 3.612739e-01, 3.687931e-01, 3.733427e-01, 3.740226e-01, 3.696330e-01, 3.649135e-01, 3.624738e-01, 3.534647e-01, 3.445755e-01, 3.354665e-01, 3.239376e-01, 3.160584e-01, 3.033997e-01, 2.923008e-01, 2.808819e-01, 2.681732e-01, 2.575442e-01, 2.443856e-01, 2.379162e-01, 2.248675e-01, 2.137186e-01, 2.024198e-01];
S_unslotted_time_estimation=[0,8.189464e-02, 1.546945e-01, 2.192381e-01, 2.801020e-01, 3.279672e-01, 3.724228e-01, 4.081392e-01, 4.418358e-01, 4.703630e-01, 4.953505e-01, 5.145885e-01, 5.256874e-01, 5.393361e-01, 5.475552e-01, 5.497250e-01, 5.536546e-01, 5.526847e-01, 5.537146e-01, 5.522648e-01, 5.473853e-01, 5.380062e-01, 5.312869e-01, 5.234177e-01, 5.145085e-01, 5.032697e-01, 4.929907e-01, 4.838516e-01, 4.709529e-01, 4.594241e-01, 4.442756e-01];
%S_unslotted_time_estimation_channel_estimation=[0,8.349000e-02, 1.502800e-01, 2.061100e-01, 2.517000e-01, 2.876500e-01, 3.125900e-01, 3.330800e-01, 3.477200e-01, 3.546600e-01, 3.580400e-01, 3.609300e-01, 3.565900e-01, 3.513900e-01, 3.431600e-01, 3.353200e-01, 3.235100e-01, 3.121100e-01, 3.034300e-01, 2.897700e-01, 2.768200e-01, 2.639700e-01, 2.507700e-01, 2.367900e-01, 2.293200e-01, 2.126200e-01, 2.035000e-01, 1.910200e-01, 1.804500e-01, 1.687700e-01, 1.591000e-01];
%S_unslotted_time_estimation_channel_estimation=[0, 7.941206e-02, 1.437656e-01, 1.958004e-01, 2.347565e-01, 2.635436e-01, 2.837016e-01, 3.016998e-01, 3.082992e-01, 3.120488e-01, 3.123788e-01, 3.085691e-01, 3.031797e-01, 2.958404e-01, 2.877712e-01, 2.782022e-01, 2.681132e-01, 2.565143e-01, 2.432457e-01, 2.300570e-01, 2.201380e-01, 2.069193e-01, 1.958004e-01, 1.830717e-01, 1.734727e-01, 1.633137e-01, 1.555144e-01, 1.427757e-01, 1.327367e-01, 1.261174e-01, 1.166183e-01];
S_unslotted_time_estimation_channel_estimation=[0,7.953205e-02, 1.497450e-01, 2.054295e-01, 2.546845e-01, 2.953705e-01, 3.286171e-01, 3.561644e-01, 3.766223e-01, 3.913709e-01, 4.046795e-01, 4.084192e-01, 4.093191e-01, 4.119788e-01, 4.107189e-01, 4.033297e-01, 4.026097e-01, 3.936106e-01, 3.826317e-01, 3.749625e-01, 3.654135e-01, 3.531447e-01, 3.463554e-01, 3.325767e-01, 3.155284e-01, 3.081992e-01, 2.954305e-01, 2.813319e-01, 2.705129e-01, 2.580742e-01, 2.481752e-01];
S_slotted_ALOHA_ideal=G.*exp(-G);
S_unslotted_ALOHA_ideal=G.*exp(-2.*G);

Time_error=[2.250840e-02 4.660194e-02 6.812825e-02 9.115695e-02 1.136853e-01 1.316498e-01 1.560812e-01 1.749844e-01 1.969231e-01 2.187440e-01];

figure;
plot(G,S_slotted_ALOHA, '--+', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [1 0 0]); hold on;
%plot(G, S_slotted_ALOHA_GDMA, '-+', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [1 0 0]);
plot(G, S_unslotted_ALOHA, '--o', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 1]);hold on;
%plot(G, S_unslotted_ALOHA_GDMA, '-o', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 1]);
plot(G,S_slotted_ALOHA_ideal, ':', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 0]);
plot(G,S_unslotted_ALOHA_ideal, ':', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 0]);
%plot(G,S_unslotted_time_estimation,':o', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 1 1])

set(gca, 'fontsize', 13);
xlabel('$G$','interpreter','latex', 'fontsize', 23);
ylabel('Throughput', 'fontsize', 21);
%leg = legend('Slotted_ALOHA_PerfectCSI', 'Slotted_ALOHA_GDMA_PerfectCSI','UnSlotted_ALOHA_PerfectCSI_PerfectTime', 'UnSlotted_ALOHA_GDMA_PerfectCSI_PerfectTime','UnSlotted_ALOHA_GDMA_PerfectCSI_Estimate_Time','Slotted_ALOHA_Theory','Unslotted_ALOHA_Theory' ,'location', 'southwest');
leg = legend('Slotted_ALOHA_PerfectCSI', 'UnSlotted_ALOHA_PerfectCSI_PerfectTime', 'Slotted_ALOHA_Theory','Unslotted_ALOHA_Theory' ,'location', 'southwest');
set(leg, 'Interpreter', 'latex');
set(leg, 'fontsize', 12);
set(gca, 'fontname', 'Times New Roman');
set(gcf,'position',[100, 100, 850, 600]);
grid on;


figure;
plot(G,S_slotted_ALOHA_GDMA_es, '-+', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [1 0 0]); hold on;
plot(G, S_slotted_ALOHA_GDMA, '-+', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 1]);

set(gca, 'fontsize', 13);
xlabel('$G$','interpreter','latex', 'fontsize', 23);
ylabel('Throughput', 'fontsize', 21);
leg = legend('Slotted_ALOHA_GDMA_BlindlyEstimation', 'Slotted_ALOHA_GDMA_PerfectCSI','location', 'southeast');
set(leg, 'Interpreter', 'latex');
set(leg, 'fontsize', 12);
set(gca, 'fontname', 'Times New Roman');
set(gcf,'position',[100, 100, 850, 600]);
grid on;

figure;
%plot(G, S_unslotted_ALOHA, '--o', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 1]);hold on;
plot(G, S_unslotted_ALOHA_GDMA, '-o', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 1]);hold on;
plot(G,S_unslotted_time_estimation,'-d', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [1 0.5 0])
plot(G,S_unslotted_time_estimation_channel_estimation,'-s', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [1 0 0])
plot(G,S_unslotted_ALOHA_ideal, ':', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 0]);
%plot(G,S_slotted_ALOHA_ideal, '--', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 0]);


set(gca, 'fontsize', 13);
xlabel('$G$','interpreter','latex', 'fontsize', 23);
ylabel('Throughput', 'fontsize', 21);
leg = legend( 'UnSlotted_ALOHA_GDMA_PerfectCSI_PerfectTime','UnSlotted_ALOHA_GDMA_EstimatedChannel_PerfectTime','UnSlotted_ALOHA_GDMA_EstimatedChannel_EstimatedTime','Unslotted_ALOHA_Theory' ,'Slotted_ALOHA_Theory' ,'location', 'southwest');
set(leg, 'Interpreter', 'latex');
set(leg, 'fontsize', 12);
set(gca, 'fontname', 'Times New Roman');
set(gcf,'position',[100, 100, 850, 600]);
grid on;

%{
figure;
plot([0.1:0.1:1],Time_error, '--o', 'MarkerSize', 12, 'linewidth', 1.5, 'Color', [0 0 1]);
set(gca, 'fontsize', 13);
xlabel('$G$','interpreter','latex', 'fontsize', 23);
ylabel('Time error', 'fontsize', 21);
set(leg, 'Interpreter', 'latex');
set(leg, 'fontsize', 12);
set(gca, 'fontname', 'Times New Roman');
set(gcf,'position',[100, 100, 850, 600]);
grid on;
%}
