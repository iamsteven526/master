clc 
clear

c_x=dlmread('c_x.txt');
c_y=dlmread('c_y.txt');
cc_x=dlmread('cc_x.txt');
cc_y=dlmread('cc_y.txt');
figure

scatter(c_x,c_y);
hold on;
scatter(cc_x,cc_y)