% Dispersion scan retrieval algorithm
% Based on "Fast iterative retrieval algorithm for ultrashort pulse
% characterization using dispersion scans" paper
clc
clear

%% Load and process data
load('wavelength.mat') %spectrometer wavelength calibration file

addpath('G:\Atto\Data\LASC\MHz\mhzharmonics\2019-02-25\');% file with scans on the atto server


img = imread('siscan1ststage3','jpg');
% figure(1);
% imshow(img);
% img(:,1:225) = [];
% img(img<5) = 0;

figure(1);
for k = 1:2
img = imdiffusefilt(img); % some anisotropic filtering to the image
end
z = linspace(-1,3,300); % mm glass insertion range. think about proper calibration
n = 4;
wl = arrayfun(@(i) mean(wl(i:i+n-1)),1:n:length(wl)-n+1); %resizing the wavelength vector to match the binned image

colormap(jet)
img = imresize(imrotate(img,11,'bilinear','crop'),0.25); %4x binning, final image size 300x480
imagesc(wl,z,flipud(img)) %measured scan
set(gca,'YDir','normal')
ylabel('glass,mm')
xlabel('wavelength, nm')


%% First guess
%Constants
tau = 7e-15; %FWHM of the pulse in time domain in s
c = 3e8; %speed of light in vacuum, m/s
wl_0 = 800e-9; %central wavelength of the fundamental, m. 
w0 = 2*pi*(c/wl_0); %central frequency of the fundamental, rad/s
sigma = 8*log(2)/(tau^2); %std.dev 1/s^2


A = 1; %amplitude of gaussian guess

%spectral phase guess
GVD = 50.621; %BK7 group velocity dispersion, fs^2/mm
[maxvalue,idx] = max(img); 

GDD = GVD.*z(end-idx); %fs^2


w = 2*pi*c./(wl.*1e-9); %frequency of sh in rad/s
% plot(w,GDD);

wl_1 = 680e-9; %fundamental wavelength range
wl_2 = 920e-9;



 w2 = linspace(2*pi*c./(wl_1),2*pi*c./(wl_2),480);
% 
% GDD = interp1(w2,GDD,w2,'linear');

% plot(w2,GDD);
% 

GD = cumtrapz(w2,GDD);
phase = cumtrapz(w2,GD);

%  plot(w2,phase)
 


% %Gaussian pulse guess
Gauss = A.*exp(-((w2-w0).^2)./sigma).*exp(1i.*phase);
% 
% plot(w2,(abs(Gauss)).^2)


%% RETRIEVE THIS B

%phase matrix
wl = (2*pi*c./w2).*1e-9;
n=sqrt(1+1.03961212./(1-0.00600069867./wl.^2)+0.231792344./(1-0.0200179144./wl.^2)+1.01046945./(1-103.560653./wl.^2)); %BK7 refractive index

k = w2.*n./c;

zk = (z*1e-3)'*k;

GGuess = kron(Gauss,ones(300,1)); 
iter = 1;

pm = exp(1i.*zk);

%multiply phase matrix by your gaussian guess
% for k = 1:10


% plot(w2,(abs(Gauss)).^2)

GGuess = GGuess.*pm; %Initial guess multiplied by pm


%to time domain
% t = linspace(2*pi/w2(1),2*pi/w2(end),480);
Gt = ifft(GGuess,[],2);

%shg
Gshgt = Gt.^2;

%to frequency domain
Gshgw = fft(Gshgt,[],2);

%multiply by measured intensity
Gup = sqrt(im2double(img(:,:))).*exp(1i.*angle(Gshgw));
figure(2)
imagesc(abs(Gup))

%to time domain
Gtr = ifft(Gup,[],2);

%P coeff
P = Gtr.*conj(Gt);

%next guess
Gt2 = nthroot(abs(P),3).*exp(1i.*angle(P));

%to frequency domain and remove phase from glass

Gw2 = fft(Gt2,[],2).*conj(pm);

deltaZ = z(2)-z(1);


% GGuess = sum(Gw2)./deltaZ;
% end

imagesc(abs(Gw2))

