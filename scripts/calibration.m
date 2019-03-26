%calibration


%% read the image and convert RGBa to RGB and grayscale
img = imread('spectrometercalibration','jpg');
addpath('\\fysfile01\atomhome$\iv2837sy\My Documents\MATLAB\siscan\calibration')
load('matlab') %reference spectra for the diodes from avantes

%camera image
figure(2);
imshow(img);


%% processing
image = img';
img_slice = image(:,519:554)';
PixelSp = sum(img_slice);

% figure(3);
% plot(PixelSp);
p = polyfit([1077 1595],[ 404.8 451.1],1);
x1 = 1:1920;
y1 = polyval(p,x1);

figure(4);
plot(x1,y1);
title('Fitting Curve')
xlabel('Pixel')
ylabel('Wavelength')


figure(5);
title('Spectrum Comparison')
plot(y1,PixelSp,'b');
hold on
plot(wl,(I405-dark405)./10);
plot(wl,(I450-dark)./10);
xlabel('\lambda, nm')
ylabel('Intensity, a.u.')
legend('Dscan','405ref','450ref','Location','northwest')
xlim([315 470])
ylim([0 6000])

wavelength = 'wavelength.mat';
wl = y1;
save(wavelength,'wl');

