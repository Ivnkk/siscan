
load('wavelength.mat')

img = imread('siscan1ststage2','jpg');
figure(1);
imshow(img);
img(:,1:225) = [];

figure(2);
for k = 1:2
img = imdiffusefilt(img);
end
insert = linspace(2,-2,800);
colormap(jet)
imagesc(wl,insert,img)
set(gca,'YDir','normal')
ylabel('glass,mm')
xlabel('wavelength, nm')

save