function scan = processScan(img,fund,z,wl)
%PROCESSSCAN processes taken image from the camera into dscan that can be
%accepted by the algorithm
% requires fundamental spectrum
% glass insertion - from file

%   Detailed explanation goes here


% fundamental parameters
wf = fund(:,1)';
c = 3e8;
%spectrometer calibrated wavelength binning
n = 4;
wl = arrayfun(@(i) mean(wl(i:i+n-1)),1:n:length(wl)-n+1);
w_exp = 2*pi*c./(wl.*1e-9)*1e-15; % to frequency (second harmonic from calib)

img = imresize(img,0.25,'nearest');% resizing the image

w3 = wf+wf(end); %second harmonic from fundamental
img = double(img);

[W_exp,Z]= meshgrid(w_exp,z); %dublicate the vectors
[W3,Z2]=meshgrid(w3,z);

img = interp2(W_exp,Z,img,W3,Z2,'spline',0);

scan = img;

end

