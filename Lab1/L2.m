clear variables;
close all;
clc;

z = tf('z');

H = 0.3/(z-0.7)

nyquist(H)
bode(H)

figure;
imp = zeros(10,1);
imp(1) = 1;
B = 0.3; A = [1, -0.7];
oi = filter(B,A,imp);
freqz(B,A);
