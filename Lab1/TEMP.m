clear variables;
close all;
clc;

x = [1,0,0,0,1,0,0,0];

fft_x = fft(x);

gwm_x = (abs(fft(x)).^2) / 8;

plot((0:7), gwm_x);

