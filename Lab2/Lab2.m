%% PDF z zadaniami

clear variables;
close all;
clc;

z = tf('z');
z_1 = 1/z;

AA = [1 -0.7];
AB = 1;
BA = [0 0.3];
BB = [1 -1];

HA = 0.3 * z_1 / (1-0.7 * z_1);
HB = 1 - z_1;

figure;
hold on;
nyquist(HA);
nyquist(HB);
hold off;

figure;
freqz(BA,AA);
figure
freqz(BB,AB);


N = 0:1:9;
[~, N_max] = size(N);
imp = zeros(N_max,1);
imp(1) = 1;

yA = filter(BA,AA,imp);
yB = filter(BB,AB,imp);

figure;
hold on;
stairs(N, yA);
stairs(N, yB);
hold off;

% Biały szum
N = 2000;
e = randn(N, 1);

syg1 = filter(BA, AA, e);
syg2 = filter(BB, AB, e);

figure;
plot(1:N, syg1);
figure;
plot(1:N, syg2);

GWM_we = 20*log10(abs(fft(e,N)/N));
GWM1_wy = 20*log10(abs(fft(syg1,N)/N));
GWM2_wy = 20*log10(abs(fft(syg2,N)/N));

figure;
hold on;
plot(1:N, GWM_we);
plot(1:N, GWM1_wy);
hold off;

figure;
hold on;
plot(1:N, GWM_we);
plot(1:N, GWM2_wy);
hold off;

%% Laborki
clear variables;
close all;
clc; 

% -------------------------------------------------------------------------
% Zadanie 3.1 - Generacja sygnalu
% -------------------------------------------------------------------------

N = 16000; % Liczba probek sygnalu dyskretnego
i = 0:1:N;

A = 0.25; % Amplituda sinusow
f_1 = 200; % Czestotliwosc pierwszego sinusa [Hz]
f_2 = 400; % Czestotliwosc drugiego sinusa [Hz]
f_3 = 800; % Czestotliwosc trzeciego sinusa [Hz]

f_s = 8000; % Czestotliwosc probkowania [Hz]

y_1 = A * sin(2 * pi * f_1/f_s * i); % Sygnal sinusoidalny pierwszy
y_2 = A * sin(2 * pi * f_2/f_s * i); % Sygnal sinusoidalny drugi
y_3 = A * sin(2 * pi * f_3/f_s * i); % Sygnal sinusoidalny trzeci

x_1 = y_1 + y_2 + y_3; % Suma sygnalow sinusoidalnych

% Wykreslenie sygnalu x_1 (sumy sinusow) i obciecie osi X
figure;
plot(i / f_s, x_1);
xlim([0 2*(1/min([f_1 f_2 f_3]))]);

% GWM z rozdzielczoscia 1 Hz
figure;
fft_x = fft(x_1,f_s)/N;
gwm_x = 20*log10(abs(fft_x));
plot((0:f_s/2),gwm_x(1:f_s/2+1));
xlabel('f [Hz]')

% -------------------------------------------------------------------------
% Zadanie 3.2 - Decymacja
% -------------------------------------------------------------------------

% Wybor co K-tej probki z sygnalu x_1
K = 4; % Wybieramy co K-ta probke
x_2 = x_1(1:K:N+1);
[~, K_max] = size(x_2);

figure;
plot(0:N/K, x_2);

fft_1 = fft(x_1(N-(200-1):N), f_s)/N;
gwm_1 = 20*log10(abs(fft_1));
fft_2 = fft(x_2(K_max-(200-1):K_max), f_s)/K_max;
gwm_2 = 20*log10(abs(fft_2));

figure;
hold on;
plot((0:f_s/2),gwm_1(1:f_s/2+1)); % Do cz. Nyquista
plot((0:f_s/2),gwm_2(1:f_s/2+1)); % Do cz. Nyquista
hold off;
xlabel('f [Hz]')