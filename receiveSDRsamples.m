%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSVB
% SDR DAB Project
% sample RF signal
% HuT, HSLU-T&A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

%% 1. parameter definition

fc = 227.36e6;        % center frequency [Hz]
rs = 1.792e6;         % sampling rate [Hz]
framelen = 256*20;    % frame length [samples]
nframes = 1000;       % number of frames to be received
Tcp = 246*10^-6;      % Länge des cycle Peräfix
Tu = 1000*10^-6;      % Länge eines Symbols
N = rs*Tu;            % number of Sampels for baseband

% create receiver object
hSDRrRx = comm.SDRRTLReceiver(...
    'CenterFrequency', fc, ...
    'EnableTunerAGC',  true, ...
    'SampleRate',      rs, ...
    'SamplesPerFrame', framelen, ...
    'OutputDataType',  'double');

% sample
R = zeros(framelen,nframes);
for i = 1:nframes
    [y,len,lost,late] = step(hSDRrRx);
    if (lost>0)
        error('error: samples lost');
    end
    R(:,i) = y;
end

display spectrum
hSpectrum = dsp.SpectrumAnalyzer(...
    'Name',             'Power Density Spectrum',...
    'Title',            'Power Density Spectrum', ...
    'SpectrumType',     'Power density',...
    'FrequencySpan',    'Full', ...
    'SampleRate',       rs, ...
    'YLimits',          [-130,0],...
    'SpectralAverages', 50, ...
    'FrequencySpan',    'Start and stop frequencies', ...
    'StartFrequency',   -rs/4, ...
    'StopFrequency',    rs/4,...
    'Position',         figposition([50 30 30 40]));
for i = 1:nframes
    step(hSpectrum, R(:,i));
end

r = R(:);
a_coff1 = getPhaseRefVector(1);
A = [zeros(ceil((N-length(a_coff1))/2),1); a_coff1 ; zeros(floor((N-length(a_coff1))/2),1)];

figure(1)
scatter(imag(A),real(A))
grid on;
title("Pole/Zero plot of the coefficients")
xlabel("Real (A)")
ylabel("Imag (A)")
figure(2)
plot(-length(A)/2:1:length(A)/2-1,abs(A))
title("Transmission function from the coefficients")
xlabel("Frequency [Hz]")
ylabel("Abs Value")


% Transform the Vector into the time Domain
 a = ifft(fftshift(A));

% Correlate the Signal to get the Peaks
a_corr = abs(conv(r(end:-1:1)',a));
%a_corr = filter(s(end:-1:1),1,r);
a_corr = a_corr(end:-1:1);
figure(3)
plot(1:length(a_corr),a_corr)
title("Correlation of vector a and r")
xlabel("Sample")
ylabel("Value")

% Find the Position of the peaks
frame_length = round(length(a_corr)/30);
for ii=1:30
[~,pos(ii)] = max(a_corr((ii-1)*frame_length+1:(ii)*frame_length));
peaks(ii) = pos(ii)+(ii-1)*frame_length;
end    
figure(4)
plot(1:length(a_corr),a_corr)
hold on
plot(peaks,a_corr(peaks),'*')
title("Correlation of vector a and r, with Peak detection")
xlabel("Sample")
ylabel("Value")
legend('Corr. Signal','Peaks','Location','south')


%% 4 Chanel estimation
position = peaks(10)-N+1;
Y_10 = fft((r(position:position+N-1)));
H_10 = fftshift(Y_10).*conj(A);
h_10 = ifft(fftshift(H_10));
figure(5)
plot(abs(h_10))
title("Impulse response")
xlabel("Sample")
ylabel("Value")
figure(8)   
plot(abs(H_10))
title("Frequency response")
xlabel("Frequency [Hz]")
ylabel("Value")
Y_11 = fft(r(peaks(11):peaks(11)+N-1))
H_11 = fftshift(Y_11).*conj(A);
h_11 = ifft(fftshift(H_11));
figure(6)
plot(abs(h_11))
hold on
plot(abs(h_10-h_11))



