% method1: ?????
clear all;
close all;
clc;

[s1, Fs] = audioread('1200.wav');
[s2, Fs] = audioread('1500.wav');
[s3, Fs] = audioread('1800.wav');
[s4, Fs] = audioread('2100.wav');
[s5, Fs] = audioread('2400.wav');

len = [length(s1), length(s2), length(s3), length(s4), length(s5)];
samples_min = min(len);
s1 = s1(1:samples_min);
s2 = s2(1:samples_min);
s3 = s3(1:samples_min);
s4 = s4(1:samples_min);
s5 = s5(1:samples_min);

in = input('Please enter the current speed of the car:','s');
speed = str2num(in);
sam_d = mod(speed,5); 

%% synthesis
if speed <= 5 && speed > 0
    gate = 0;
    signal = s2;
    samples = floor(samples_min/((1200 + sam_d*60)/1500));
elseif speed > 5 && speed <= 10
    gate = 0;
    signal = s3;
    samples = floor(samples_min/((1500 + sam_d*60)/1800));
elseif  speed > 10 && speed <= 15
    gate = 0;
    signal = s4;
    samples = floor(samples_min/((1800 + sam_d*60)/2100));
elseif  speed > 15 && speed <= 20
    gate = 0;
    signal = s5;
    samples = floor(samples_min/((2100 + sam_d*60)/2400));
elseif  speed == 0 % idlt
    gate = 0;
    signal = s1;
else 
    fprintf('error');
    gate = 1;
end

if gate == 0
    if sam_d==0
        s_out = signal; 
    else 
        xx = 1:1:samples_min;
        xx2 = linspace(1, samples_min, samples);
        s_out = interp1(xx,signal,xx2,'spline');
        s_out = s_out(1:samples_min);
    end

    % Fade in and fade out
    alpha = 0.99;
    a = ((alpha+1)/(1-alpha))^0.5;
    endstep = 48000;
    step = -endstep:endstep;
    k = log(a)*22050/(endstep);
    Fadein = 0.5*(tanh(k*(step-1)/22050)+1)';
    Fadeout = 1 - Fadein;
    Filter = [Fadein;ones(samples_min-length(Fadein)-length(Fadeout), 1);Fadeout];
    Filter = transpose(Filter);
    s_out = s_out.*Filter;

    sound(s_out,Fs)
    audiowrite('output.wav', s_out, Fs); 

%%
% spectrum and time domain plot
    x=signal;
    fx=fft(x);
    figure;
    subplot(211);
    plot(x);
    subplot(212)
    N=length(x);
    f=(0:N-1)*Fs/N;
    % disg=20*log10(abs(fx));
    % plot(f(1:N/2+1),disg(1:N/2+1))
    plot(f(1:N/2+1),abs(fx(1:N/2+1)));

    x=s_out;
    fx=fft(x);
    figure;
    subplot(211);
    plot(x);
    subplot(212)
    N=length(x);
    f=(0:N-1)*Fs/N;
    % disg=20*log10(abs(fx));
    % plot(f(1:N/2+1),disg(1:N/2+1))
    plot(f(1:N/2+1),abs(fx(1:N/2+1)));
end
% % spectrum
% figure;       
% win_sz = 128;       
% han_win = hanning(win_sz);      
% nfft = 2^nextpow2(length(han_win));       
% nooverlap = win_sz - 1;
% [S, F, T, P] = spectrogram(s2, window, nooverlap, nfft, Fs);   
% imagesc(T, F, log10(abs(S)))  
% colorbar;
% set(gca, 'YDir', 'normal') 
% xlabel('Time (secs)')
% ylabel('Freq (Hz)')
% title('STFT transform spectrum')