%Implement in Matlab a stereo three-band tone control with switchable loudness control.

clear all
close all

Fs = 44100; %sampling frequency
input = [1;zeros(1000,1)]; %impulse signal
lowCutoffFrequency = 250 / (Fs/2) * pi; %convert to radian frequency 
highCutoffFrequency = 4000 / (Fs/2) * pi; %convert to radian frequency

for toneLowValue = 0:0.25:1 %sweeps through knob values of low-shelving filter
    for toneMidValue = 0:0.25:1 %sweeps through knob values of mid-shelving filter
        for toneHighValue = 0:0.25:1 %sweeps through knob values of high-shelving filter
            output = toneControl(input, lowCutoffFrequency, highCutoffFrequency, toneLowValue, toneMidValue, toneHighValue, 0);
            B = output; %B = toneControl() output
            [X, w] = freqz(B); %frequency response (X), angular frequency (w), for the output of tone Control
            F = w / pi * Fs/2; %creates a frequency vector
            fontsize = 20;

            xrange = [20 20000];
            subplot(2,1,1)
            semilogx(F, 20*log10(abs(X))); %plots the magnitude response of
            title(sprintf('lowpot: %f, midpot: %f, highpot: %f', toneLowValue, toneMidValue, toneHighValue)) %formats the title to show changing values for low, mid and high knobs.
            g = gca(); %sets g to current axis, subplot(2,1,1)
            g.Title.FontSize = fontsize; %sets title fontsize
            g.FontSize = fontsize; %sets all other fontsizes
            xlabel('Frequency[Hz] (Magnitude Response)')
            ylabel('Gain[dB]')
            xlim(xrange); %sets limits for xaxis

            subplot(2,1,2)
            semilogx(F, angle(X)); %plots the phase response of frequency response F
            g = gca();
            g.FontSize = fontsize;
            xlabel('Frequency[Hz] (Phase Response)')
            ylabel('Phase [Deg]')
            xlim(xrange);

            print('-dpng', sprintf('output_for_highpot_%f_midpot_%f.png', toneHighValue, toneMidValue));
        end
    end
end

%%

%Test code for audio playback

[y, Fs] = audioread('longer-sample.wav');

lowCutoffFrequency = 1000/(Fs/2) * pi;
highCutoffFrequency =  5000/(Fs/2) * pi;
toneLowValue = 0.5;
toneMidValue = 0.5;
toneHighValue = 0.5;
loudness = 1;

% % sweep the mid pot control
for n = 0:1
    toneLowValue = n * 0.25; 
    output = toneControl(y, lowCutoffFrequency, highCutoffFrequency, toneLowValue, toneMidValue, toneHighValue, loudness);

    plot(y)
    player = audioplayer(y, Fs);
    player.play();
    pause
    
    plot(output)
    player = audioplayer(output, Fs);
    player.play();
    pause
    
end

%%
function output = toneControl(input, lowCutoffFrequency, highCutoffFrequency, toneLowValue, toneMidValue, toneHighValue, loudnessValue)

%assume I have already converted to radian frequency when passing values to this function

%lowCutoffFrequency, highCutoffFrequency

B = highCutoffFrequency - lowCutoffFrequency;

wcLs = lowCutoffFrequency;
wcHs = highCutoffFrequency;
wcPn = sqrt(wcLs * wcHs);

%translate knob values to gain values 

LsG = 10^((toneLowValue * 24 - 12)/20);

PnG = 10^((toneMidValue * 24 - 12)/20);

HsG = 10^((toneHighValue * 24 - 12)/20);

%%
%low shelving filter

tLs = tan(wcLs/2); %creates a variable for tan(wc/2)
s1 = sqrt(LsG); %creates a variable for sqrt(G)

zrsLs = [ (LsG * tLs + s1), (LsG * tLs - s1) ];

plsLs = [ (tLs + s1), (tLs - s1) ];
%%

%peak/notch filter

tPn = tan(B/2);
cPn = cos(wcPn); %in this case wc is the center frequency, as the p/n has two cutoffs, low and high
sPn = sqrt(PnG);

zrsPn = [ (sPn + PnG * tPn), (-2 * sPn * cPn), (sPn - PnG * tPn) ];
 
plsPn = [ (sPn + tPn), (-2 * sPn * cPn), (sPn - tPn) ];

%%

%high shelving filter

sHs = sqrt(HsG);
tHs = tan(wcHs/2);

zrsHs = [ (sHs * tHs + HsG), (sHs * tHs - HsG) ];

plsHs = [ (sHs * tHs + 1), (sHs * tHs - 1)];

%%

%filters in series

lsOut = filter(zrsLs, plsLs, input);
pnOut = filter(zrsPn , plsPn, lsOut);
hsOut = filter(zrsHs, plsHs, pnOut);

load('a_weighting_FIR_coefficients','bn');

% loudness control

if loudnessValue == 1
    
output = filter(bn, 1, hsOut);
    
elseif loudnessValue == 0
    
output = hsOut;

else

fprints('expected Boolean value')

end

end

