clear variables
close all

x = 1:50;
data = x + rand(1, 50)*10;

coarsening_length = 5;
r = 10;
yfilted = Gaussianfilter(data, r, coarsening_length);

plot(x, data, x, yfilted)