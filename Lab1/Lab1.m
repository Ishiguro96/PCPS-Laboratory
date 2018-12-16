close all;
clc
 
f = 400;
fp = 8000;
N = fp;
A = 0.24;
i = (0:N-1)';
x = A * sin(2 * pi * f / fp * i);
plot(i(1:100)/fp,x(1:100));
xlabel('t [s]')
 
 
figure;
% obliczenie gÄ™stoÅ›ci widmowej mocy -- rozdzielczoÅ›Ä‡ widmowa 1Hz
fft_x = fft(x,fp)/N;
gwm_x = 20*log10(abs(fft_x));
plot((0:fp/2),gwm_x(1:fp/2+1));
xlabel('f [Hz]')
 
 
figure;
 
% kwantyzacja
%b = 2;
%b = 4;
b = 4;
q = 2/(2^b);
xq = round(x/q)*q;
 
% kwantyzacja z ditheringiem
%d = 0.5*(rand(N,1)-rand(N,1));
%xq = round(x/q+d)*q;
stem((1:100),xq(1:100));
 
 
% obliczenie gÄ™stoÅ›ci widmowej mocy -- rozdzielczoÅ›Ä‡ widmowa 1Hz
figure;
fft_x = fft(xq,fp)/N;
gwm_x = 20*log10(abs(fft_x));
plot((0:fp/2),gwm_x(1:fp/2+1), 'o');
xlabel('f [Hz]')
