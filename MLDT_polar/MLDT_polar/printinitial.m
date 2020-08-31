
clear all
x = dlmread('rx_x_intitial.txt');
y = dlmread('rx_y_intitial.txt');
centroid_x= dlmread('centroid_x_intitial.txt');
centroid_y= dlmread('centroid_y_intitial.txt');
temp_x=dlmread('chCoef_x_intitial.txt');
temp_y=dlmread('chCoef_y_intitial.txt');
%t_x=dlmread('temp_x.txt');
%t_y=dlmread('temp_y.txt');


figure;
plot(x,y,'k.');
hold on;
scatter(centroid_x,centroid_y,'rx')

hold on;
scatter(temp_x,temp_y,'o')
%scatter(t_x,t_y,'<')
legend('receive symbol','estimation','chCoef','Location','southeast');
grid on;