% Dispersion scan retrieval algorithm
% Based on "Fast iterative retrieval algorithm for ultrashort pulse
% characterization using dispersion scans" paper. equation numbers follow
% that of the paper.
%% Initialization
clc
clear

%Constants
tau = 7; %FWHM of the pulse in time domain in fs

c = 3e8; %speed of light in vacuum, m/s

wl_0 = 800e-9; %central fundamental wavelength, m. 

w0 = 2*pi*(c/wl_0)*1e-15; %central frequency, rad/fs

sigma = 8*log(2)/(tau^2); %std.dev in gaussian pulse 1/fs^2

wf_1 = 3.1473; %border of the freq. range #1

wf_2 = 1.2; %border of the freq. range #2

w2 = linspace(wf_1,wf_2,480); %frequency range of the fundamental

GVD = 50.621; %BK7 group velocity dispersion, fs^2/mm

z = linspace(-4,4,300); % glass insertion range.

iter = 1; % number of iterations in the algorithm

wl_fun = (2*pi*c./w2).*1e-9; %convert freq. to wl to calculate refr.index 

n = sqrt(1+1.03961212./(1-0.00600069867./wl_fun.^2)+0.231792344./(1-0.0200179144./wl_fun.^2)+1.01046945./(1-103.560653./wl_fun.^2)); %BK7 refractive index

k = w2.*1e15.*n./c; %wavenumber

zk = (z*1e-3)'*k; %kz matrix

N = 7; %number of iterations in the algorithm

Err = 1; %initial value for error

%% Load and process the data

%load
addpath('G:\Atto\Data\LASC\MHz\mhzharmonics\2019-02-25\');

load('wavelength.mat') %spectrometer wavelength calibration file

% load('fundamental.mat')

img = imread('siscan1ststage3','jpg'); %scan to load

%binning
n = 4; %bin factor

wl = arrayfun(@(i) mean(wl(i:i+n-1)),1:n:length(wl)-n+1); %resizing the wavelength vector to match the binned image

w_exp = 2*pi*c./(wl.*1e-9)*1e-15; %wavelength to frequency for calibrated spectrometer

img = imresize(imrotate(img,10,'bilinear','crop'),0.25); %4x binning, final image size 300x480

%reinterpolate the scan onto frequency axis
figure(1);
colormap(jet)

w3 = linspace(w_exp(1),w_exp(end),480); 
img = double(img);

[W_exp,Z]=meshgrid(w_exp,z); %dublicate the vectors
[W3,Z2]=meshgrid(w3,z);

img = interp2(W_exp,Z,img,W3,Z2,'linear');
imagesc(w3,z,img)

ylabel('glass,mm')
xlabel('frequency, rad/fs')
%% First guess

%spectral phase
[maxvalue,idx] = max(img); 

GDD = GVD.*z(end-idx); 

GDD = interp1(w2,GDD,w2,'linear');

GD = cumtrapz(w2,GDD);

phase = cumtrapz(w2,GD);

% Gaussian pulse 
Gauss = exp(-((w2-w0).^2)./sigma).*exp(1i.*phase); % eq(6)



%% RETRIEVE

GGuess = kron(Gauss,ones(300,1)); %extend gauss pulse vector into matrix

pm = exp(1i.*zk); %phase matrix. eq(7)

for i = 1:N
GGuess = GGuess.*pm; %Initial guess multiplied by pm. eq(8)

%frequency shift (if needed)
dw = w2(end)-w2(end-1);
N = 480;
t = linspace(0,(N-1)*2*pi/dw/N,480); 

Gt = ifft(GGuess,[],2); %to time domain. eq(9)

Gshgt = Gt.^2; %generate shg. eq(10)

Gshgw = fft(Gshgt,[],2); %to frequency domain. eq(11)

Gup = sqrt(double(img)).*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)

Gtr = ifft(Gup,[],2); %to time domain. eq(13)

P = Gtr.*conj(Gt); %P coeff. eq(14)

Gt2 = nthroot(abs(P),3).*exp(1i.*angle(P)); %next guess. eq(15)

Gw2 = fft(Gt2,[],2).*conj(pm); %to frequency domain and remove phase from glass. eq(16) and (17)

GGuess = sum(Gw2.*(z(2)-z(1)))./(z(end)-z(1)); %new field. eq(18)

%error estimation

retr = abs(Gshgw).^2./max(max(abs(Gshgw).^2)); %normalized retrieved scan

img1 = img./max(max(img));

mu = sum(img1.*retr)./sum(retr); %minimization vector/constant

Err = sum(sum((img1 - mu.*retr).^2))./(length(z)*length(w2));

%plot the trace
figure(2);
colormap(jet)
imagesc(w_exp,z,abs(Gshgw).^2)
ylabel('glass,mm')
xlabel('frequency, rad/fs')
drawnow;

end

%% PLOT THE RESULTS

figure(3);
plot(wl_fun,abs(GGuess).^2./max(abs(GGuess).^2.))
xlim([0.6 1])
ylabel('Spectral power, a.u.')
xlabel('\lambda, nm')

