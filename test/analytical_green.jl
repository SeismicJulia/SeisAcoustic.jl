function [d_time]=GH_greens_funct(rd,dt,nt,c,src,opt)

% Inputs:
% f: frequency [Hz]
% rd: zeros(ns,nr);src and rev distance [m].
% src: zeros(1,nt);usually Ricker.
%
% Outputs:
% d_time: data (Greens function) in the time domain.
%
%  Copyright (C) 2019, Exploration and Environmental Geophysics (EEG) at
%  ETH Zurich

src=reshape(src,[1,nt]);

ff = 1/dt;    % [Hz] Sampling rate
fs = linspace(0,ff,nt); % [Hz] frequency vector

k = 2*pi*fs/c;
if strcmp(opt,'greens1d')
    greens = @(r,f) -1i./(2*k) .* exp(-1i*r*k);
elseif strcmp(opt,'greens2d_1st')
    greens = @(r,f) 1i/4 * besselh(0,2, r*k) .* (-2*pi*fs*1i); % first-order system
elseif strcmp(opt,'greens2d_2nd')
    greens = @(r,f) 1i/4 * besselh(0,2, r*k); % second-order system
elseif strcmp(opt,'greens3d')
    greens = @(r,f) exp(-1i*r*k)./(4*pi*r);
else
    error('opt not supported.\n');
end

%% ---------
SRC=fft(src);

d_f = zeros(length(fs),length(rd));
for l=1:length(rd)
    d_f(:,l) = SRC .* sum( greens(rd(:,l),fs) , 1 );
end
d_f( isnan(d_f) ) = 0;

%% 6. Display in time domain
d_time = ifft( d_f, [], 1, 'symmetric' );


end
