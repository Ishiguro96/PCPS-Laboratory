%{
    LABORATORIUM PCPS
    DAWID TOBOR
%}

%% ZADANIE 1
clear variables;
close all;
clc;

fp = 500; % Czestotliwosc probkowania
f = 40; % Czestotliwosc sygnalu
A = 1; % Amplituda sygnalu 
N = fp;

i = (0:1000) / fp; % Wektor dyskretny czasu

% Prostokat
figure;
wykresProstokat = A * square(2 * pi * f * i);
stairs(i, wykresProstokat);
title('Przebieg czasowy - prostokat');
ylabel('x');
xlabel('t[s]');

% Sinus
figure;
wykresSinus = A * sin(2 * pi * f * i);
stairs(i, wykresSinus);
title('Przebieg czasowy - sinus');
ylabel('x');
xlabel('t[s]');

%Kwantyzacja
kwant = uencode(wykresSinus, 2);
kwant = udecode(kwant, 2);
figure;
plot(i, kwant);

% Bialy szum
figure;
wykresBialySzum = randn(1001, 1);
plot(i, wykresBialySzum);
title('Przebieg czasowy - biały szum');
ylabel('x');
xlabel('t[s]');

fft_x = fft(wykresSinus, fp)/N;
gwm_x = 20*log10(abs(fft_x));
figure;
plot((0:fp/2), gwm_x(1:fp/2+1));

%% zad 2
x = importdata('../Materiały/syg_2018/sygnaly_00/syg06_3000Hz');
fp = 7000;
N = fp;
fn=fp/2;
pr=[0:1:fp-1];
X = fft(x,fp);
GWM = (abs(X).^2)/N
plot(0:fn,10*log10(GWM(1:fn+1)))
title('GWM')
ylabel('dB');
xlabel('f[Hz]');

fft_x = fft(x, fp)/N;
gwm_x = 20*log10(abs(fft_x));
figure;
plot((0:fp/2), gwm_x(1:fp/2+1));

figure;
fp=10000;
t=(0:99)/fp;
plot(t,x(1:100));
xlabel('t[s]');
title('przebieg czasowy');
ylabel('Amplituda');

%% zad 3
N=fp;
fn=fp/2;
pr=[0:1:fp-1];
%GWM dla sin
x=sin(2*pi*fx/fp*pr);
X=fft(x,fp);
GWM=(abs(X).^2)/N
plot(0:fn,10*log10(GWM(1:fn+1)))
title('GWM dla sin')
ylabel('dB');
xlabel('f[Hz]');

%GWM dla dla prost
x=square(2*pi*fx/fp*pr)
X=fft(x,fp);
GWM=(abs(X).^2)/N
title('GWM dla prost')
ylabel('dB');
xlabel('f[Hz]');
plot(0:fn,(GWM(1:fn+1))) %dla prostokatnego w skali liniowej

%GWM dla bia�ego sumu
x=randn(fp,1)
X=fft(x,fp);
GWM=(abs(X).^2)/N
plot(0:fn,10*log10(GWM(1:fn+1)))
title('GWM dla szumu')
ylabel('dB');
xlabel('f[Hz]');

%% zad 4 
fx = 50;
N = fp;
fp=10000;
pr=[0:1:fp-1];
x=sin(2*pi*fx/fp*pr);
t=(0:fp-1)/fp;
 
kwant = uencode(x,8);
kwant = udecode(kwant,8);
plot(t,kwant);
figure;
X=fft(kwant,fp);
GWM=(abs(X).^2)/N
plot(0:fx,10*log10(GWM(1:fx+1)))

%% zad 5

fp = 100000; % cz�stotliwo�� pr�bkowania (Hz)
fn = fp / 2;
L = 125; % liczba pr�bek
t = (0:L-1)/fp; % Wektor czasu

x5 = sin(2*pi*250*t) + sin(2*pi*450*t) + sin(2*pi*1000*t) + sin(2*pi*2033*t) + sin(2*pi*4900*t);
x51 = decimate(x5, 2);
x52 = decimate(x5, 3);
x53 = decimate(x5, 4);

plot((fp*t), x51)
title('Wykres sygna�u wielosinusoidalnego')
xlabel('czas dyskretny [pr�bki]')
ylabel('warto�� sygna�u')
grid on

X = fft(x51, fp);
GWM = abs(X).^2 / L;
plot(0:fn, 10*log10(GWM(1:fn+1)))
title('GWM')
xlabel('f [Hz]')
ylabel('dB')
grid on


