function scan = processScan(img,fund,z,wl)
%PROCESSSCAN processes taken image from the camera into dscan that can be
%accepted by the algorithm
% requires fundamental spectrum
% glass insertion - from file

%   Detailed explanation goes here


% fundamental parameters
f = fund(:,1)';

%spectrometer calibrated wavelength binning

wl = fliplr(wl);

wl = imresize(wl,1/3.75,'bilinear');
f_exp = 299.792./wl;

img = fliplr(img);

img = imresize(imresize(img,[1024 2048],'bilinear'),0.25,'bilinear');% resizing the image

f_shg = f + min(f); %second harmonic from fundamental
img = double(img);

[F_exp,Z]= meshgrid(f_exp,z); %dublicate the vectors
[F_shg,Z2]=meshgrid(f_shg,z);

img = interp2(F_exp,Z,img,F_shg,Z2,'spline',0);

scan = img;

end

