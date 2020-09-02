
clear all
x = dlmread('rx_x_intitial.txt');
y = dlmread('rx_y_intitial.txt');
centroid_x= dlmread('centroid_x_intitial.txt');
centroid_y= dlmread('centroid_y_intitial.txt');
temp_x=dlmread('chCoef_x_intitial.txt');
temp_y=dlmread('chCoef_y_intitial.txt');

tx=[-0.828177 -0.496076 -1.94693 0.505825 0.618864 1.91973 -0.575543 0.895126];
ty=[-0.699546 1.58674 0.954138 -1.55632 -0.0982055 -0.944299 0.167358 0.784126];

figure;
plot(x,y,'k.');
hold on;
scatter(centroid_x,centroid_y,'rx')

hold on;
scatter(temp_x,temp_y,'o')
%scatter(tx,ty,'<')

legend('receive symbol','estimation','chCoef','Location','southeast');
grid on;