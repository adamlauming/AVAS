clear all;
close all;
clc;

[s1, Fs] = audioread('sampleslib/1200.wav');
[s2, Fs] = audioread('sampleslib/1500.wav');
[s3, Fs] = audioread('sampleslib/1800.wav');
[s4, Fs] = audioread('sampleslib/2100.wav');
[s5, Fs] = audioread('sampleslib/2400.wav');

len = [length(s1), length(s2), length(s3), length(s4), length(s5)];
minimumsamples = min(len);
samples_min = 96000;

% fade in and fade out
alpha = 0.99;
a = ((alpha+1)/(1-alpha))^0.5;
endstep = 12000;
step = (-endstep+1):endstep;
k = log(a)*22050/(endstep);
Fadein = 0.5*(tanh(k*(step-1)/22050)+1)';
Fadeout = 1 - Fadein;
Filter = [Fadein;ones(samples_min-length(Fadein)-length(Fadeout), 1);Fadeout];
Filter = transpose(Filter);

% fade in
alpha = 0.99;
a = ((alpha+1)/(1-alpha))^0.5;
endstep = 12000;
step = -endstep:endstep;
k = log(a)*22050/(endstep);
Fadein = 0.5*(tanh(k*(step-1)/22050)+1)';
Filter_in = [Fadein;ones(samples_min-length(Fadein), 1)];
Filter_in = transpose(Filter_in);

% fade out
alpha = 0.99;
a = ((alpha+1)/(1-alpha))^0.5;
endstep = 12000;
step = -endstep:endstep;
k = log(a)*22050/(endstep);
Fadein = 0.5*(tanh(k*(step-1)/22050)+1)';
Fadeout = 1 - Fadein;
Filter_out = [ones(samples_min-length(Fadeout), 1);Fadeout];
Filter_out = transpose(Filter_out);

speedlist = xlsread('speed.xlsx');
s_output = [ones(length(speedlist)*samples_min, 1)];
s_output = transpose(s_output);
iteration = 0; % times of the same speed

for index = 1:length(speedlist)
    speed = speedlist(index)
    if speed <= 20
        sam_d = mod(speed,5); 
        if index == 1
            s_out = s1(1:samples_min);
            % Fade in
            s_out = s_out.*Filter_in;

        elseif index >= 2 && index < length(speedlist)
            speed_before = speedlist(index-1);
            s_out_before = s_out;
            if speed == speed_before
                if (iteration+1)*samples_min > minimumsamples
                    iteration = 0;
                else 
                    iteration = iteration + 1;
                end
            else
                iteration = 0;
            end

            %% synthesis
            if speed <= 5 && speed > 0
                signal = s2(iteration*samples_min+1:(iteration+1)*samples_min);
                samples = floor(samples_min/((1200 + sam_d*60)/1500));
            elseif speed > 5 && speed <= 10
                signal = s3(iteration*samples_min+1:(iteration+1)*samples_min);
                samples = floor(samples_min/((1500 + sam_d*60)/1800));
            elseif  speed > 10 && speed <= 15
                signal = s4(iteration*samples_min+1:(iteration+1)*samples_min);
                samples = floor(samples_min/((1800 + sam_d*60)/2100));
            elseif  speed > 15 && speed <= 20
                signal = s5(iteration*samples_min+1:(iteration+1)*samples_min);
                samples = floor(samples_min/((2100 + sam_d*60)/2400));
            elseif  speed == 0 % idlt ??
                signal = s1(iteration*samples_min+1:(iteration+1)*samples_min);
            end

            if sam_d==0
                s_out = signal; 
            else 
                xx = 1:1:samples_min;
                xx2 = linspace(1, samples_min, samples);
                s_out = interp1(xx,signal,xx2,'spline');
                s_out = s_out(1:samples_min);
            end

            if speed == speed_before || speed_before > 20
                s_out = signal;
            else 
                % cross fade in and fade out
                s_out = s_out.*Filter_in;
                s_out_before = s_out_before.*Filter_out;
                s_out(1:24000) = s_out(1:24000) + s_out_before(72001:96000); 
            end

        elseif index == length(speedlist)
            speed_before = speedlist(index-1);
            s_out_before = s_out;
            if speed == speed_before
                s_out = s1(1:samples_min);
                % Fade out
                s_out = s_out.*Filter_out;
            else 
                s_out_before = s_out_before.*Filter_out;
                s_out = s1(1:samples_min);
                % cross Fade in and fade out
                s_out = s_out.*Filter;
                s_out(1:24000) = s_out(1:24000) + s_out_before(72001:96000); 
            end
        end
    else % more than 20km/h
        s_out = zeros(1,samples_min);
    end 
    
    s_output(((index-1)*samples_min+1):index*samples_min) = s_out;
    audiowrite('output.wav', s_output, Fs); 
end

x=s_output;
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

%% get the sound pressure level of speech signal
x = s_output;
Length=length(x);

framlen = 100;
M=Fs*framlen/1000;
m = mod(Length,M);
if m >= M/2 
    x = [x,zeros(1,M-m)];
    Length = length(x);
else   
    l = floor(Length/M);
    x = x(1,1:M*l);
    Length = length(x);
end
N = Length/M;
s = zeros(1,N);
spl = zeros(1,N);
for k = 1:N
    s =x((k-1)*M + 1:k*M);
    spl(1,k) = SPLCal(s,Fs,framlen);
end
% plot
t = 1:Length;
SPL = zeros(1,Length);
for r = 1:N
    SPL(1,(r-1)*M+1:r*M) = spl(r);
end
figure;
subplot(211)
plot(t/Fs,x);
grid on
xlabel('Time(s)');
title('Output signal');
subplot(212)
stairs(t/Fs,SPL,'r');
grid on
xlabel('Time(s)');
ylabel('Sound pressure level(dB)');
title('Sound pressure level of speech signal(dB)/100ms');