clear;
close all
load('colorm.mat');
load('fund1.mat');
n1=480; % number of points in the frequency domain

n0=256; % number of points in the glass insertion domain

f0=linspace(0.2,0.8,n1)+eps; %frequency vector

cf=0.375; %central frequency

tau2 = 7; % pulse duration

glass=linspace(-3,3,n0); % glass insertion vector

wl=300./(f0); % frequency to wavelength

f1 = fnspec(:,1)/2/pi;

% inten=fnspec(:,2)'.*exp(-(tau2*(f0-cf)).^2); % gaussian intensity profile

inten=interp1(f1,fnspec(:,2)',f0,'linear',0);
inten=inten./max(inten); % normalize

id=f0<cf; %step?
t = linspace(-1,1,n1);
% [~,phase0]= 200*((f0-cf)/0.25).^2.*(~id); % phase of the simulated field
[~,phase0]= generate_phase(linspace(900,700,480),800,[133 690 -200 80 pi]);
% phase0 = phase0  + 0.8.*rand(1,length(f0));
phase0 = phase0.*0.02.*sin(2*pi*4*t) + 0.0*gallery('normaldata',size(f0),1);
fshg=f0 + min(f0); %shg frequency

dphase=2*pi*wextend('ac','sp0',glass',length(f0)-1,'r').*nBK7(wextend('ar','sp0',wl,length(glass)-1,'d')/1000)./(wextend('ar','sp0',wl,length(glass)-1,'d')/1000000); % phase matrix

intshg=SHGv21(dphase,inten,phase0); %create an SHG

imagesc(fshg,glass,intshg); %plot

colormap(cmap0); %cmp

% xlim([0.1,1]);
% % 
save('sample0.mat','fshg','glass','intshg','f0','inten');
