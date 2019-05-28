clear;

load('colorm.mat');

n1=1024; % number of points in the frequency domain

n0=256; % number of points in the glass insertion domain

f0=linspace(0.00000,1.0,n1)+eps; %frequency vector

glass=linspace(-6,6,n0); % glass insertion vector

wl=300./(f0); % frequency to wavelength

wl2 = linspace(wl(1),wl(end),length(wl));

inten=exp(-(15*(f0-0.375)).^2); % intensity profile

inten=inten./max(inten); % normalize

cf=0.375; %central frequency

id=f0<cf; %step?

phase0=200*((f0-cf)/0.25).^2.*(~id); % phase of the simulated field

fshg=f0+min(f0); %shg frequency

dphase=2*pi*wextend('ac','sp0',glass',length(f0)-1,'r').*nBK7(wextend('ar','sp0',wl,length(glass)-1,'d')/1000)./(wextend('ar','sp0',wl,length(glass)-1,'d')/1000000); % phase matrix

intshg=SHGv21(dphase,inten,phase0); %create an SHG

imagesc(f0,glass,intshg); %plot

colormap(cmap0); %cmp

xlim([0.5,1]);

save('sample0.mat','fshg','glass','intshg','f0','inten','phase0');
