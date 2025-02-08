clear all
close all
[x1 x2]=meshgrid(-3:0.1:3);
f = x1.^4 - 3*x1.^2 + x2.^2 - 2*x1 -2*x2 +10;
surf(x1,x2,f)
axis([-2 2 -2 2 4 20])