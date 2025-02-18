function [Ib, T, cutoff_gamma, cutoff_s] = loggaborconv2(im, filterbankargs)
% ========================================================================
% Returns the filter response results of an image "im" going through several
% log-gabor filters contained in a specific filter bank given whose arguments 
% as "filterbankargs".
% 
% Input: (1) im - the grayscale image to pass the filters.
%        (2) filterbankargs - the structure contains necessary fields to
%        shape mainly the band and orientation density of the filter bank.
%        Its fields include nscale, norient, minWaveLength, mult, sigmaOnf
%        and dThetaOnSigma.
% Output: (1) Ib - the convolution results through each of the filters
%         contained in cell array form. The index of the cell array is
%         arranged as (s+(o-1)*nscale) where o indicates in which 
%         orientation the passed filter has been limit and s indicates
%         which out of the nscale bands in the same orientation this filter
%         is located.
%         (2) T - the noise compensation factor.
%         (3) cutoff_gamma - the gain factor that controls the sharpness
%         of the cutoff of the filter response spread.
%         (4) cutoff_s - the cutoff value of the filter response spread.
%         
%
% This function is an adaptation of Peter Kovesi's PHASECONG2.
% ========================================================================
% PHASECONG2 - Computes edge and corner phase congruency in an image.
%
% Copyright (c) 1996-2009 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby  granted, free of charge, to any  person obtaining a copy
% of this software and associated  documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% The software is provided "as is", without warranty of any kind.
% References:
%
%     Peter Kovesi, "Image Features From Phase Congruency". Videre: A
%     Journal of Computer Vision Research. MIT Press. Volume 1, Number 3,
%     Summer 1999 http://mitpress.mit.edu/e-journals/Videre/001/v13.html

[nscale, norient, minWaveLength, mult, sigmaOnf, dThetaOnSigma] = checkfilterbankargs(filterbankargs);

cutoff_gamma = 1./sigmaOnf;
cutoff_s = 1./minWaveLength;

[rows,cols] = size(im);
imagefft = fft2(im);              % Fourier transform of image
thetaSigma = pi/norient/dThetaOnSigma;

zero = zeros(rows,cols);
if mod(cols,2)
    xrange = [-(cols-1)/2:(cols-1)/2]/(cols-1);
else
    xrange = [-cols/2:(cols/2-1)]/cols;
end
if mod(rows,2)
    yrange = [-(rows-1)/2:(rows-1)/2]/(rows-1);
else
    yrange = [-rows/2:(rows/2-1)]/rows;
end
[x,y] = meshgrid(xrange, yrange);
radius = sqrt(x.^2 + y.^2);
theta = atan2(-y,x);
lp = ones([rows,cols]);
rs = 0.45;
sl = 0.05;
lp(radius>rs) = exp(-(radius(radius>rs)-rs).^2/(2*sl^2));
lp = ifftshift(lp);
radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters
theta  = ifftshift(theta);        % are constructed with 0 frequency at the corners.
radius(1,1) = 1;                  % Get rid of the 0 radius value at the 0
% frequency point (now at top-left corner)
% so that taking the log of the radius will
% not cause trouble.
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;    % save a little memory


logGabor = cell(1,nscale);
Ib = cell(1,nscale);
ifftFilterArray = cell(1,nscale);

for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));
    logGabor{s} = logGabor{s}.*lp;        % Apply low-pass filter
    logGabor{s}(1,1) = 0;                 % Set the value at the 0 frequency point of the filter
    % back to zero (undo the radius fudge).
    
    ifftFilterArray{s} = real(ifft2(logGabor{s}))*sqrt(rows*cols);
end

spread = cell(1,norient);

if norient == 1
    EM_n = sum(sum(logGabor{1}.^2));
    for s = 1:nscale
        Ib{s} = 2.*real(ifft2(imagefft .* logGabor{s}));
    end
else
    for o = 1:norient
        angl = (o-1)*pi/norient;           % Filter angle.
        
        % For each point in the filter matrix calculate the angular distance from
        % the specified filter orientation.  To overcome the angular wrap-around
        % problem sine difference and cosine difference values are first computed
        % and then the atan2 function is used to determine angular distance.
        
        ds = sintheta * cos(angl) - costheta * sin(angl);    % Difference in sine.
        dc = costheta * cos(angl) + sintheta * sin(angl);    % Difference in cosine.
        dtheta = abs(atan2(ds,dc));                          % Absolute angular distance.
        spread{o} = exp((-dtheta.^2) / (2 * thetaSigma^2));  % Calculate the
        % angular filter component.
    end
    for o = 1:norient                    % For each orientation.        
        for s = 1:nscale                  % For each scale.
            filter = logGabor{s} .* spread{o};   % Multiply radial and angular
            % components to get the filter.
            
            ifftFilt = real(ifft2(filter))*sqrt(rows*cols);  % Note rescaling to match power
            ifftFilterArray{s} = ifftFilt;                   % record ifft2 of filter
            % Record mean squared filter value at smallest
            EM_n = sum(sum(filter.^2));           % scale. This is used for noise estimation.
            Ib{s+(o-1)*nscale} = 2.*real(ifft2(imagefft .* filter));
            
        end
    end
end

% Compensate for noise
% We estimate the noise power from the energy squared response at the
% smallest scale.  If the noise is Gaussian the energy squared will have a
% Chi-squared 2DOF pdf.  We calculate the median energy squared response
% as this is a robust statistic.  From this we estimate the mean.
% The estimate of noise power is obtained by dividing the mean squared
% energy value by the mean squared filter value

medianE2n = median(reshape(Ib{1}.^2,1,rows*cols));
meanE2n = -medianE2n/log(0.5);

noisePower = meanE2n/EM_n;                       % Estimate of noise power.

EstSumAn2 = zero;
for s = 1:nscale
    EstSumAn2 = EstSumAn2 + ifftFilterArray{s}.^2;
end

EstSumAiAj = zero;
for si = 1:(nscale-1)
    for sj = (si+1):nscale
        EstSumAiAj = EstSumAiAj + ifftFilterArray{si}.*ifftFilterArray{sj};
    end
end
sumEstSumAn2 = sum(sum(EstSumAn2));
sumEstSumAiAj = sum(sum(EstSumAiAj));

k = 2.0;
EstNoiseEnergy2 = 2*noisePower*sumEstSumAn2 + 4*noisePower*sumEstSumAiAj;

tau = sqrt(EstNoiseEnergy2/2);                     % Rayleigh parameter
EstNoiseEnergy = tau*sqrt(pi/2);                   % Expected value of noise energy
EstNoiseEnergySigma = sqrt( (2-pi/2)*tau^2 );

T =  EstNoiseEnergy + k*EstNoiseEnergySigma;       % Noise threshold

T = T/1.7;        % Empirical rescaling of the estimated noise effect to
% suit the PC_2 phase congruency measure

function [nscale, norient, minWaveLength, mult, sigmaOnf, dThetaOnSigma] = checkfilterbankargs(filterbankargs)
if isfield(filterbankargs, 'nscale')
    nscale = filterbankargs.nscale;
else
    nscale = 4;
end
if isfield(filterbankargs, 'norient')
    norient = filterbankargs.norient;
else
    norient = 4;
end
if isfield(filterbankargs, 'minWaveLength')
    minWaveLength = filterbankargs.minWaveLength;
else
    minWaveLength = 3;
end
if isfield(filterbankargs, 'mult')
    mult = filterbankargs.mult;
else
    mult = 2.1;
end
if isfield(filterbankargs, 'sigmaOnf')
    sigmaOnf = filterbankargs.sigmaOnf;
else
    sigmaOnf = 0.55;
end
if isfield(filterbankargs, 'dThetaOnSigma')
    dThetaOnSigma = filterbankargs.dThetaOnSigma;
else
    dThetaOnSigma = 1.5;
end