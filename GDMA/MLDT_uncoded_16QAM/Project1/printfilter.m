clear all

syn = dlmread('syn.txt');
drift=dlmread('drift.txt');
drift_real=dlmread('drift_real.txt');
drift_img=dlmread('drift_img.txt');
syn_real=dlmread('syn_real.txt');
syn_img=dlmread('syn_img.txt');
figure;
plot([0.25:0.25:20],syn,'-o');
hold on;
plot([0.25:0.25:20],drift,'-d');
legend('synchronization','time drift:0.48Tc','Location','southwest');
xlabel('chip time');
grid on;
figure;
subplot(2,1,1)
plot(drift_real,'-o');
subplot(2,1,2)
plot(syn_real,'-d');

figure;
subplot(2,1,1);
plot(drift_i);