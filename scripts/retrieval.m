% Dispersion scan retrieval algorithm
% Based on "Fast iterative retrieval algorithm for ultrashort pulse
% characterization using dispersion scans" paper. equation numbers follow
% that of the paper.
clc
clear

%% Load the data 

addpath('G:\Atto\Data\LASC\MHz\mhzharmonics\2019-02-25\');

load('wavelength.mat') % imaging spectrometer wavelength calibration file

load('fundamental.mat') % fundamental spectrum

Int = abs(Int)./max(Int); %normalize fundamental

img = imread('siscan1ststage3','jpg'); %scan to load


%% Constants

tau = 6; %FWHM of the pulse in time domain in fs

c = 3e8; %speed of light in vacuum, m/s

wl_0 = 800e-9; %central fundamental wavelength, m. 

w0 = 2*pi*(c/wl_0)*1e-15; %central frequency, rad/fs

sigma = 8*log(2)/(tau^2); %std.dev in gaussian pulse 1/fs^2

wf_1 = 2*pi*(c/wlf(1))*1e-6 +0.4; %border of the freq. range #1

wf_2 = (2*pi*(c/wlf(end))*1e-6); %border of the freq. range #2

wf = linspace(wf_1,wf_2,480); %frequency range of the fundamental

wf3 = (2*pi*(c./wlf)*1e-6);

Int = interp1(wf3,Int,wf,'linear',0)./wf.^2;

Int = abs(Int)./max(Int);

GVD = 50.621; %BK7 group velocity dispersion, fs^2/mm

z = linspace(-4,4,300); % glass insertion range.

iter = 1; % number of iterations in the algorithm

wl_fun = (2*pi*c./wf).*1e-9; %convert freq. to wl to calculate refr.index 

n = sqrt(1+1.03961212./(1-0.00600069867./wl_fun.^2)+0.231792344./(1-0.0200179144./wl_fun.^2)+1.01046945./(1-103.560653./wl_fun.^2)); %BK7 refractive index

k = wf.*1e15.*n./c; %wavenumber

zk = (z*1e-3)'*k; %kz matrix

N = 18; %number of iterations in the algorithm

Err = 1; %initial value for error

%% Process the data



%binning
n = 4; %bin factor

wl = arrayfun(@(i) mean(wl(i:i+n-1)),1:n:length(wl)-n+1); %resizing the wavelength vector to match the binned image

w_exp = 2*pi*c./(wl.*1e-9)*1e-15; %wavelength to frequency for calibrated spectrometer

img = imresize(imrotate(img,11,'bilinear','crop'),0.25); %4x binning, final image size 300x480

%reinterpolate the scan onto frequency axis
figure(1);
colormap(jet)
w3 = wf+wf(end);
% w3 = linspace(w_exp(1),w_exp(end),length(w_exp)); 
img = double(img);

[W_exp,Z]= meshgrid(w_exp,z); %dublicate the vectors
[W3,Z2]=meshgrid(w3,z);

img = interp2(W_exp,Z,img,W3,Z2,'linear',0);
imagesc(w3,z,img)

ylabel('glass,mm')
xlabel('frequency, rad/fs')
%% First guess

%spectral phase
[maxvalue,idx] = max(img); 

GDD = GVD.*z(end-idx); 

GDD = interp1(wf,GDD,wf,'linear');

GD = cumtrapz(wf,GDD);

phase = cumtrapz(wf,GD) + wf.*110;

% Gaussian pulse 
Gauss = exp(-((wf-w0).^2)./sigma).*exp(1i.*phase); % eq(6)



%% RETRIEVE

GGuess = kron(Gauss,ones(300,1)); %extend gauss pulse vector into matrix

pm = exp(1i.*zk); %phase matrix. eq(7)

Meas = sqrt(double(img)); %measured d scan field

for i = 1:N
    


GGuess = GGuess.*pm; %Initial guess multiplied by pm. eq(8)

%frequency shift (if needed)
dw = abs(mean(diff(wf)));
Np = length(w_exp);
t = linspace(0,(Np-1)*2*pi/dw/Np,Np); 

Gt = ifft(GGuess,[],2); %to time domain. eq(9)

Gshgt = Gt.^2; %generate shg. eq(10)

Gshgw = fft(Gshgt,[],2); %to frequency domain. eq(11)

Gup = Meas.*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)

Gtr = ifft(Gup,[],2); %to time domain. eq(13)

P = Gtr.*conj(Gt); %P coeff. eq(14)

Gt2 = nthroot(abs(P),3).*exp(1i.*angle(P)); %next guess. eq(15)

Gw2 = fft(Gt2,[],2).*conj(pm); %to frequency domain and remove phase from glass. eq(16) and (17)

GGuess = sum(Gw2.*(z(2)-z(1)))./(z(end)-z(1)); %new field. eq(18)

% GGuess = sqrt(Int).*exp(1i.*angle(GGuess)); %multiply by fundamental

%error estimation

retr = abs(Gshgw).^2./max(max(abs(Gshgw).^2)); %normalized retrieved scan

img1 = img./max(max(img)); %normalized measured

% img1 = fillmissing(img1,'constant',0);

mu = sum(img1.*retr)./sum(retr)+eps; %minimization vector/constant

% mu = fillmissing(mu,'constant',0);

Err = sum(sum((img1 - mu.*retr).^2))./(length(z)*length(wf)); %error estimation

% Meas = Meas./mu;

%plot the trace
figure(2);
colormap(jet)
imagesc(w3,z,abs(Gshgw).^2)
ylabel('glass,mm')
xlabel('frequency, rad/fs')
drawnow;

end

%% PLOT THE RESULTS
% GGuess = GGuess.*exp(-1i.*wf*0000);
phase = angle(GGuess);
phase = unwrap(phase) ;

%spectral domain
figure(3);
plot(wf,abs(GGuess).^2./max(abs(GGuess).^2.))


yyaxis left
ylabel('Spectral power, a.u.')
xlabel('freq, rad/fs')
yyaxis right
plot(wf,phase)
ylabel('Spectral phase, rad')


%temporal domain
% dt = 2*pi*(wf(1)-wf(end))./wf.^2;

% t = linspace(-dt*(Np-1),dt*(Np-1),Np);



figure(4);
% [t,Ret] = IFFT(wf,GGuess);
% Ret = ifftshift(Ret);
Ret = fft(GGuess);
phase = unwrap(angle(Ret));


plot(t,abs(Ret).^2./max(abs(Ret).^2))

yyaxis left
ylabel('Intensity, a.u.')
xlabel('time, fs')
yyaxis right
plot(t,phase)
ylabel('Temporal phase, rad')







