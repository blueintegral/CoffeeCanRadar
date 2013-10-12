function userinput

sprintf('How would you like to process the data?')
prompt = '1 for RealTime Audio Processing 2 for First acquiring data and then processing';
result = input(prompt);

if result == 1
    final_workingrealtime;
userinput;
elseif result == 2
    notrealtime;
userinput;
else
    sprintf('Undefined user input')
end


%MIT IAP Radar Course 2011
%Resource: Build a Small Radar System Capable of Sensing Range, Doppler, 
%and Synthetic Aperture Radar Imaging 
%
%Gregory L. Charvat

%Process Range vs. Time Intensity (RTI) plot

%NOTE: set up-ramp sweep from 2-3.2V to stay within ISM band
%change fstart and fstop bellow when in ISM band
function final_workingrealtime
clear all;
close all;
pause on
%read the raw data .wave file here
Yt =  audiorecorder(16000,16,2);
% Y = wavrecord(44100*5,44100,2);
FS = 16000;
NBITS = 16;

%constants
c = 3E8; %(m/s) speed of light

%radar parameters
Tp = 20E-3; %(s) pulse time
N = Tp*FS; %# of samples per pulse
fstart = 2260E6; %(Hz) LFM start frequency for example
fstop = 2590E6; %(Hz) LFM stop frequency for example
%fstart = 2402E6; %(Hz) LFM start frequency for ISM band
%fstop = 2495E6; %(Hz) LFM stop frequency for ISM band
BW = fstop-fstart; %(Hz) transmti bandwidth
f = linspace(fstart, fstop, N/2); %instantaneous transmit frequency
x = 0;
%range resolution
rr = c/(2*BW);
max_range = rr*N/2;
record(Yt);
pause(8)
n = 1;


%%%%
%Find the size of 



%%%
while(1)  
Ybar = getaudiodata(Yt, 'single');
sz = size(Ybar);
%the input appears to be inverted
%Method 2
if x == 0
    Y = Ybar;
    x = 1;
else
    Y = [Y(size(Ybar,1)-5000:size(Y,1),:);Ybar(5000:size(Ybar,1),:)];
end
trig = -1*Y(:,1);
s = -1*Y(:,2);



%parse the data here by triggering off rising edge of sync pulse
count = 0;
thresh = 0;
start = (trig > thresh);
for ii = 100:(size(start,1)-N)
    if start(ii) == 1 & mean(start(ii-11:ii-1)) == 0
        %start2(ii) = 1;
        count = count + 1;
        sif(count,:) = s(ii:ii+N-1);
        time(count) = ii*1/FS;
    end
end
%subtract the average
ave = mean(sif,1);
for ii = 1:size(sif,1);
    sif(ii,:) = sif(ii,:) - ave;
end

zpad = 8*N/2;

figure(20);
sif2 = sif(2:size(sif,1),:)-sif(1:size(sif,1)-1,:);
v = ifft(sif2,zpad,2);
S=v;
R = linspace(0,max_range,zpad);
for ii = 1:size(S,1)
    %S(ii,:) = S(ii,:).*R.^(3/2); %Optional: magnitude scale to range
end
S = dbv(S(:,1:size(v,2)/2));
m = max(max(S));
imagesc(R,time,S-m,[-80, 0]);
colorbar;
ylabel('time (s)');
xlabel('range (m)');
title('RTI with 2-pulse cancelor clutter rejection');


record(Yt)
pause(2)
stop(Yt)
end
end




function notrealtime
%read the raw data .wave file here
Y =  wavrecord(16000*10,16000,2);
% Y = wavrecord(44100*5,44100,2);
FS = 16000;
NBITS = 16;

%constants
c = 3E8; %(m/s) speed of light
%radar parameters
Tp = 20E-3; %(s) pulse time
N = Tp*FS; %# of samples per pulse
fstart = 2260E6; %(Hz) LFM start frequency for example
fstop = 2590E6; %(Hz) LFM stop frequency for example
%fstart = 2402E6; %(Hz) LFM start frequency for ISM band
%fstop = 2495E6; %(Hz) LFM stop frequency for ISM band
BW = fstop-fstart; %(Hz) transmti bandwidth
f = linspace(fstart, fstop, N/2); %instantaneous transmit frequency
x = 0;
%range resolution
rr = c/(2*BW);
max_range = rr*N/2;
trig = -1*Y(:,1);
s = -1*Y(:,2);
%parse the data here by triggering off rising edge of sync pulse
count = 0;
thresh = 0;
start = (trig > thresh);
for ii = 100:(size(start,1)-N)
    if start(ii) == 1 & mean(start(ii-11:ii-1)) == 0
        %start2(ii) = 1;
        count = count + 1;
        sif(count,:) = s(ii:ii+N-1);
        time(count) = ii*1/FS;
    end
end
%subtract the average
ave = mean(sif,1);
for ii = 1:size(sif,1);
    sif(ii,:) = sif(ii,:) - ave;
end

zpad = 8*N/2;

figure(20);
sif2 = sif(2:size(sif,1),:)-sif(1:size(sif,1)-1,:);
v = ifft(sif2,zpad,2);
S=v;
R = linspace(0,max_range,zpad);
for ii = 1:size(S,1)
    %S(ii,:) = S(ii,:).*R.^(3/2); %Optional: magnitude scale to range
end
S = dbv(S(:,1:size(v,2)/2));
m = max(max(S));
imagesc(R,time,S-m,[-80, 0]);
colorbar;
ylabel('time (s)');
xlabel('range (m)');
title('RTI with 2-pulse cancelor clutter rejection');
end
end