clear;
close all
load('colorm.mat');

n1=1024; % number of points in the frequency domain

n0=256; % number of points in the glass insertion domain

f0=linspace(0.25,1,n1)+eps; %frequency vector

cf=0.375; %central frequency

tau2 = 16; % pulse duration

glass=linspace(-6,6,n0); % glass insertion vector

wl=300./(f0); % frequency to wavelength

inten=exp(-(tau2*(f0-cf)).^2); % gaussian intensity profile

inten=inten./max(inten); % normalize

id=f0<cf; %step?

% [~,phase0]= 200*((f0-cf)/0.25).^2.*(~id); % phase of the simulated field
[~,phase0]= generate_phase(linspace(900,700,1024),800,[122000 600000 -620 80 0]);
fshg=f0 + min(f0); %shg frequency

dphase=2*pi*wextend('ac','sp0',glass',length(f0)-1,'r').*nBK7(wextend('ar','sp0',wl,length(glass)-1,'d')/1000)./(wextend('ar','sp0',wl,length(glass)-1,'d')/1000000); % phase matrix

intshg=SHGv21(dphase,inten,phase0); %create an SHG

imagesc(fshg,glass,intshg); %plot

colormap(cmap0); %cmp

xlim([0.1,1]);
% % 
save('sample0.mat','fshg','glass','intshg','f0','inten');
