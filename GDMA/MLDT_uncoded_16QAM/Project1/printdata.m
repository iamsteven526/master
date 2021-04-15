

clear all
x = dlmread('rx_x.txt')
y = dlmread('rx_y.txt')

ce_x= dlmread('ce_x.txt')
ce_y= dlmread('ce_y.txt')

ch_x= dlmread('ch_x.txt')
ch_y= dlmread('ch_y.txt')

centroid_x= dlmread('centroid_x.txt')
centroid_y= dlmread('centroid_y.txt')
figure;
%c = linspace(1,10,length(x));
scatter(x,y,15,'k.');

hold on;
scatter(ce_x,ce_y,50,'d')

hold on;
scatter(centroid_x,centroid_y,50,'rx')

hold on;
scatter(ch_x,ch_y,50,'o')


grid on;

legend('receive symbol','estimation','centroid','chCoef','Location','southeast');
%legend('receive symbol','chCoef','centroid','Location','southeast');

