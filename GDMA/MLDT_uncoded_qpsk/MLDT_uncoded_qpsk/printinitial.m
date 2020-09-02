clear all
x = dlmread('rx_x_intitial.txt');
y = dlmread('rx_y_intitial.txt');

centroid_x= dlmread('centroid_x_intitial.txt');
centroid_y= dlmread('centroid_y_intitial.txt');

figure;
plot(x,y,'k.');
hold on;
plot(centroid_x,centroid_y,'rx')

grid on;