% Dispersion scan retrieval algorithm
% Based on "Fast iterative retrieval algorithm for ultrashort pulse
% characterization using dispersion scans" paper


%% Load and process data
load('wavelength.mat') %spectrometer wavelength calibration file

addpath('G:\Atto\Data\LASC\MHz\mhzharmonics\2019-02-25\');% file with scans on the atto server


img = imread('siscan1ststage3','jpg');
% figure(1);
% imshow(img);
% img(:,1:225) = [];
% img(img<5) = 0;

figure();
for k = 1:2
img = imdiffusefilt(img); % some anisotropic filtering to the image
end
z = linspace(2,6,300); % glass insertion range. think about proper calibration
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
tau = 6; %FWHM of the pulse in time domain in fs
c = 3e8; %speed of light in vacuum, m/s
wl_0 = 400e-9; %central wavelength, m. 
w0 = 2*pi*(c/wl_0)*1e-15; %central frequency, rad/fs
sigma = 8*log(2)/(tau^2); %std.dev 1/fs^2

A = 1; %amplitude of gaussian guess

%spectral phase guess
GVD = 50.621; %BK7 group velocity dispersion, fs^2/mm
[maxvalue,idx] = max(img); 
figure();
GDD = GVD.*z(end-idx); %reinterpolate this


w = 2*pi*c./(wl.*1e-9)*1e-15; %frequency in rad/fs
plot(w,GDD);


 w2 = linspace(w(80),w(420),280);
% 
GDD = interp1(w,GDD,w2,'linear');
plot(w2,GDD);
% % 
GD = cumtrapz(w2,GDD);
phase = cumtrapz(w2,GD);
% % % 
% for k = 1:length(w)-1
% %     GD(k+1) = GD(k) + h*GDD
% % end

 plot(w2,phase)
% 
% 
%  
% 
% %Gaussian pulse guess
Gauss = A.*exp(-((w2-w0).^2)./sigma).*exp(1i.*phase);
% 
plot(w2,(abs(Gauss)).^2)


%% RETRIEVE THIS B

%phase matrix
wl = (2*pi*c./w2).*1e-9;
n=sqrt(1+1.03961212./(1-0.00600069867./wl.^2)+0.231792344./(1-0.0200179144./wl.^2)+1.01046945./(1-103.560653./wl.^2)); %BK7 refractive index

k = w2.*n./c;

zk = (z*1e3)'*k;

pm = exp(1i.*zk);


GGuess = Gauss.*pm; %Initial guess multiplied by pm



