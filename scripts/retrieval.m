% Dispersion scan retrieval algorithm
% Based on "Fast iterative retrieval algorithm for ultrashort pulse
% characterization using dispersion scans" paper. equation numbers follow
% that of the paper.
clc
clear
close all

%% Load the data 

addpath('G:\Atto\Data\LASC\MHz\mhzharmonics\2019-05-22\'); %rewrite paths when implementing into GUI
addpath('D:\PhD\Programmes\MATLAB\siscan\scripts\functions');
addpath('D:\PhD\Programmes\MATLAB\siscan\scripts\data');

load('wavelength.mat') % imaging spectrometer wavelength calibration file

load('fund2stage.mat') % fundamental spectrum

Int = fnspec(:,2)'; 

wf = fnspec(:,1)';

[maxval,idx] = max(Int);

img = imread('2ndstage3','png'); %scan to load

img = img - 0;
img = max(0,img);

%% Constants

tau = 10; %FWHM of the first guess in time domain in fs

c = 3e8; %speed of light in vacuum, m/s

w0 = wf(idx); %central frequency, rad/fs

w3 = wf+wf(end); %shg spectrum

sigma = 8*log(2)/(tau^2); %std.dev in gaussian pulse 1/fs^2

GVD = 50.621; %BK7 group velocity dispersion, fs^2/mm

z = linspace(-1.7,1.8,300); % glass insertion range.

wl_fun = (2*pi*c./wf).*1e-9; %convert freq. to wl to calculate refr.index 

% wl_fun = linspace(wl_fun(1),wl_fun(end),length(wl_fun));

n = sqrt(1+1.03961212./(1-0.00600069867./wl_fun.^2)+0.231792344./(1-0.0200179144./wl_fun.^2)+1.01046945./(1-103.560653./wl_fun.^2)); %BK7 refractive index

k = wf.*1e15.*n./c; %wavenumber

zk = (z*1e-3)'*k; %kz matrix

N = 400; % max number of iterations in the algorithm

Err = 1; %initial value for error

iter = 1; % counter

%% Process the data



%binning
n = 4; %bin factor

wl = arrayfun(@(i) mean(wl(i:i+n-1)),1:n:length(wl)-n+1); %resizing the wavelength vector to match the binned image

w_exp = 2*pi*c./(wl.*1e-9)*1e-15; %wavelength to frequency for calibrated spectrometer

img = imresize(img,0.25,'nearest'); %4x binning, final image size 300x480
for k = 1:2
    img = imdiffusefilt(img);
end
%reinterpolate the scan onto frequency axis
figure(1);
colormap(parula)
img = double(img);

[W_exp,Z]= meshgrid(w_exp,z); %dublicate the vectors
[W3,Z2]=meshgrid(w3,z);

img = interp2(W_exp,Z,img,W3,Z2,'spline',0);
imagesc(w3,z,img)
img = max(0,img);
img = img./max(max(img));
set(gca,'YDir','normal')

ylabel('glass,mm')
xlabel('frequency, rad/fs')
%% First guess

%spectral phase
[maxvalue,idx] = max(img); 

GDD = GVD.*z(end-idx); 

GDD = interp1(wf,GDD,wf,'linear');

GD = cumtrapz(wf,GDD);

phase = cumtrapz(wf,GD) - 50;



int_z=sum(img,2);

int_z=int_z/sum(int_z);

intzz=wextend('ac','sp0',int_z,length(wf)-1,'d');

idz=int_z==max(int_z);

Meas = sqrt(double(img)); %measured d scan field

pm = exp(1i.*zk); %phase matrix. eq(7)


% GGuess=abs(ifft(sqrt(fft(Meas(idz,:),[],2)),[],2));

% Gaussian pulse 
GGuess = exp(-((wf-w0).^2)./sigma).*exp(1i.*phase); % eq(6)



%% RETRIEVE

while Err>1e-4
    
if iter >= N
    break
end

GGuess = kron(GGuess,ones(length(z),1)); %extend pulse vector into matrix

Gt = ifft(GGuess.*pm,[],2); %to time domain. eq(8) and eq(9)

Gshgw = fft(Gt.^2,[],2); %generate shg and go to frequency domain. eq(11)

Gup = Meas.*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)

Gtr = ifft(Gup,[],2); %to time domain. eq(13)

P = Gtr.*conj(Gt); %P coeff. eq(14)

Gt2 = abs(P.^(1/3)).*exp(1i.*angle(P)); %next guess. eq(15)

Gw2 = fft(Gt2,[],2).*exp(-1i.*zk); %to frequency domain and remove phase from glass. eq(16) and (17)

GGuess = sum(Gw2.*intzz,1);

% GGuess = sum(Gw2.*(z(2)-z(1)))./(z(end)-z(1)); %new field. eq(18)

GGuess = sqrt(Int).*exp(1i.*angle(GGuess)); %multiply by fundamental

%error estimation

retr = abs(Gshgw).^2./max(max(abs(Gshgw).^2)); %normalized retrieved scan

mu = sum(img.*retr)./(sum(retr)); %minimization vector/constant

mu = kron(mu,ones(length(z),1));

Err = sum(sum((img - mu.*retr).^2))./(length(z)*length(wf)); %error estimation



iter = iter+1;

%plot the trace

figure(2);
subplot(2,1,1)
imagesc(w3,z,retr)
set(gca,'YDir','normal')
subplot(2,1,2)
imagesc(w3,z,img)
set(gca,'YDir','normal')
title(['iter=',num2str(iter), 'Error=',num2str(Err)])
ylabel('glass,mm')
xlabel('frequency, rad/fs')
drawnow;



end

%% PLOT THE RESULTS

%spectral domain
phase = angle(GGuess);
phase = unwrap(phase);

figure(3);
plot(wf,abs(GGuess).^2./max(abs(GGuess).^2.))
title('Spectral domain')
hold on
plot(wf,Int,'k--')
yyaxis left
ylabel('Spectral power, a.u.')
xlabel('freq, rad/fs')
legend('Retrieved', 'Measured')
yyaxis right

plot(wf,phase)
ylabel('Spectral phase, rad')
legend('Phase')
xlim([1.4 3.1])
hold off

%temporal domain
Np = length(wf);
dt = mean(2*pi*(wf(1)-wf(end))./wf.^2);
t = linspace(-dt*(Np-1),dt*(Np-1),Np);

figure(4);
Ret = ifft(GGuess);
Ret = ifftshift(Ret);
phase = unwrap(angle(Ret));
t2 = linspace(t(1),t(end),1920);
Inten = abs(Ret).^2./max(abs(Ret).^2);
Inten = interp1(t,Inten,t2,'pchip');
phase = interp1(t,phase,t2,'pchip');
plot(t2,Inten)
yyaxis left
title ('Temporal intensity profile')
ylabel('Intensity, a.u.')
xlabel('time, fs')
ylim([0 1.2])
yyaxis right
plot(t2,phase)
ylabel('Temporal phase, rad')
xlim([-150 150])







