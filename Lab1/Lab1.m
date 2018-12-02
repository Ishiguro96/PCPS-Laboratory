%{
    LABORATORIUM PCPS
    DAWID TOBOR
%}

clear all;
close all;
clc;

% ZADANIE 1
fp=2000; % Czestotliwosc probkowania
fx=40; % Czestotliwosc sygnalu
A=1; 
t=(0:1000)/fp; % Czas trwania sygnalu

%sygnal prostokatny
wykres_prostokat=A*square(2*pi*fx*t);
stairs(t,wykres_prostokat);
title('przebieg czasowy sygnal prostok�tny');
ylabel('x');
xlabel('t[s]');

%sygnal sinus
wykres_sinus=A*sin(2*pi*fx*t);
plot(t,wykres_sinus);
title('przebieg czasowy sygnal sin');
ylabel('x');
xlabel('t[s]');

%bialy szum
wykres_szum=randn(size(t));
plot(t,wykres_szum)
title('przebieg czasowy bia�y szum');
ylabel('x');
xlabel('t[s]');

%% zad 2
x=importdata('C:\Users\Tomasz M\Desktop\sygnaly_001\syg08_7kHz');
fp=7000;
N=fs;
fn=fp/2;
pr=[0:1:fp-1];
X=fft(x,fp);
GWM=(abs(X).^2)/N
plot(0:fn,10*log10(GWM(1:fn+1)))
title('GWM')
ylabel('dB');
xlabel('f[Hz]');

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
fp=10000;
pr=[0:1:fp-1];
x=sin(2*pi*fx/fp*pr);
t=(0:6999)/fp;
 
kwant = uencode(x,8);
kwant = udecode(kwant,8);
plot(t,kwant);
X=fft(kwant,fp);
GWM=(abs(X).^2)/N
plot(0:fn,10*log10(GWM(1:fn+1)))

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


