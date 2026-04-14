function spec = PowerSpectrum_Wavelet(lfp,nCycles,dt,frange,varargin)
% modefied from buzsaki code
% [specslope,spec] = bz_PowerSpectrumSlope(lfp,winsize,dt) calculates the
%slope of the power spectrum, a metric of cortical state and E/I balance
%see Gao, Peterson, Voytek 2016;  Waston, Ding, Buzsaki 2017
%
%INPUTS
%   lfp         
%   nCycles     cycles of wavelets
%   dt          deault 0.01
%
%   (optional)
%       'ints'      Intervals in which to calculate PSS. Default:[0 Inf]
%       'showfig'   true/false - show a summary figure of the results
%                   (default:false)
%       'saveMat'   put your basePath here to save/load
%                   baseName.PowerSpectrumSlope.lfp.mat  (default: false)
%       'IRASA'     (default: true) use IRASA method to median-smooth power
%                   spectrum before fitting
%                   (Muthukumaraswamy and Liley, NeuroImage 2018)
%       'nfreqs'    number of frequency values used to fitting (default: 200)
%                   
%
%OUTPUTS
%   specgram
%       .data       complex-valued spectrogram
%       .timestamps
%       .amp        log10-transformed amplitude of the spectrogram
%
%
%DLevenstein 2018
%with IRASA code modified from R. Hardstone & W. Munoz - 2019
%%
p = inputParser;
addParameter(p,'showfig',false,@islogical)
addParameter(p,'Redetect',false)
addParameter(p,'IRASA',false)
addParameter(p,'nfreqs',200)
addParameter(p,'ints',[0 inf])
addParameter(p,'space','log')
parse(p,varargin{:})
SHOWFIG = p.Results.showfig;
IRASA = p.Results.IRASA;
nfreqs = p.Results.nfreqs;
ints = p.Results.ints;
space = p.Results.space;

%% Calcluate spectrogram
if IRASA
    maxRescaleFactor = 2.9; %as per Muthukumaraswamy and Liley, NeuroImage 2018
    %Add the padding for edge effects of frequency smoothings (IRASA)
    padding = floor((nfreqs./2).*log10(maxRescaleFactor^2)./log10(frange(2)./frange(1)));
    actualRescaleFactor = 10.^((log10(frange(2)./frange(1)).*padding)./(nfreqs-1)); %Account for rounding
    nfreqs = nfreqs+2*padding;
    frange = frange .* [1/actualRescaleFactor actualRescaleFactor];
end

% calculate spectrogram
downsampleout = round(lfp.samplingRate.*dt);
dt = downsampleout/lfp.samplingRate; %actual dt another option is to do movmean on the amplitude...
[spec] = bz_WaveSpec(lfp,'frange',frange,'nfreqs',nfreqs,'ncyc',nCycles,...
    'downsampleout',downsampleout,'intervals',ints,'space',space);

spec.amp = log10(abs(spec.data));

%% IRASA before fitting
if IRASA
    %Frequencies inside padding 
    validFreqInds = padding+1 : nfreqs-padding;
    
    %Median smooth the spectrum
    for i_freq = validFreqInds
        inds = [i_freq-padding:i_freq-1 i_freq+1:i_freq+padding];
        resampledData(:,i_freq) = nanmedian(spec.amp(:,inds),2);
    end
    
    %Calculate the residual from the smoothed spectrum
    %spec.osci = (spec.amp(:,validFreqInds))-(resampledData(:,validFreqInds));
    power4fit = resampledData(:,validFreqInds);
    
    %Return the specgram to requested frequencies
    spec.freqs = spec.freqs(validFreqInds);
    spec.amp = spec.amp(:,validFreqInds);
    spec.data = spec.data(:,validFreqInds);
    
    spec.IRASAsmooth = power4fit;
    
    clear resampledData
else
   power4fit = spec.amp;
end

%% Fit the slope of the power spectrogram
rsq = zeros(size(spec.timestamps));
s = zeros(length(spec.timestamps),2);
yresid = zeros(length(spec.timestamps),length(spec.freqs));
for tt = 1:length(spec.timestamps)
    %Fit the line
    x = log10(spec.freqs);  y=power4fit(tt,:);
    s(tt,:) = polyfit(x,y,1);
    %Calculate the residuals (from the full PS)
    yfit =  s(tt,1) * x + s(tt,2);
    yresid(tt,:) = spec.amp(tt,:) - yfit; %residual between "raw" PS, not IRASA-smoothed
    %Calculate the rsquared value
    SSresid = sum(yresid(tt,:).^2);
    SStotal = (length(y)-1) * var(y);
    rsq(tt) = 1 - SSresid/SStotal;
end


%% Figure
if SHOWFIG
    
   %hist(specslope.data,10)
   specmean.all = mean(spec.amp,1);

   h=figure;
   h.Position = [100 100 1200 900];
   subplot(2,1,1)
   plot(log2(spec.freqs),specmean.all,'k','linewidth',2)
   LogScale('x',2)
   axis tight
   box off
   title('Periodogram')
   xlabel('Frequency')
   ylabel('log10(Power)')

   subplot(2,1,2)
   imagesc(spec.timestamps,log2(spec.freqs),spec.amp')
   title('Spectrogram')
   SpecColorRange( spec.amp );
   LogScale('y',2)
   ylabel('f (Hz)')
   xlabel('time (sec)')
   axis xy
end


end

