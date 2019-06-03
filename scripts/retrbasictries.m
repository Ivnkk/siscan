
%RETRFUN Summary of this function goes here
% retrieval function 

% inputs:
% scan - processed image of measured dscan already calibrated on frequency
% axis. 
% wf - fundamental frequency vector
% z - glass insertion range
% N - maximum number of iterations in the algorithm
% fund - measured fundamental spectrum

%outputs:
% retr  - retrieved dscan
% field - retrieved field
clc
clear
close all

%% Load the data 

addpath('G:\Atto\Data\LASC\MHz\mhzharmonics\2019-05-22\'); 
addpath('D:\PhD\Programmes\MATLAB\siscan\scripts\functions');
addpath('D:\PhD\Programmes\MATLAB\siscan\scripts\data');

load('wavelength.mat') % imaging spectrometer wavelength calibration file

load('fund1.mat') % fundamental spectrum

Int = fnspec(:,2)'./max(fnspec(:,2)'); 

wf = fnspec(:,1)';

[maxval,idx] = max(Int);

img = imread('try1','png'); %scan to load

img = img - 2;
img = max(0,img);



%% Initialization & constants

% constants
N = 200;

iter = 1;

Err = 1;

c = 3e8;

w3 = wf+min(wf); %shg frequencies

wl_f = 2*pi*300./wf; %wavelength vector. spacing?

z = linspace(-2,2,300); % glass insertion range.

phase = 2*pi*wextend('ac','sp0',z',length(wf)-1,'r').*nBK7(wextend('ar','sp0',wl_f,length(z)-1,'d')/1000)./(wextend('ar','sp0',wl_f,length(z)-1,'d')/1000000); % kz matrix 

n = 4; %bin factor

wl = arrayfun(@(i) mean(wl(i:i+n-1)),1:n:length(wl)-n+1); %resizing the wavelength vector to match the binned image

w_exp = 2*pi*c./(wl.*1e-9)*1e-15; %wavelength to frequency for calibrated spectrometer

img = imresize(img,0.25,'nearest'); %4x binning, final image size 300x480
% for k = 1:6
%     img = imdiffusefilt(img);
% end
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

Meas=sqrt(double(img)); %amplitude of the trace

int_z=sum(img,2); 

int_z=int_z/sum(int_z); %normalize

intzz=wextend('ac','sp0',int_z,length(wf)-1,'d');

idz=int_z==max(int_z);

GGuess=abs(ifft(sqrt(fft(Meas(idz,:),[],2)),[],2));
%% retrieval

while Err>0.02
    
    GGuess = kron(GGuess,ones(length(z),1)); %extend pulse vector into matrix
    
    Gt = ifft(GGuess.*exp(1i.*phase),[],2); %to time domain. eq(8) and eq(9)
    
    Gshgw = fft(Gt.^2,[],2); %generate shg and go to frequency domain. eq(11)
    
    Gup = Meas.*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)
    
    Gtr = ifft(Gup,[],2); %to time domain. eq(13)
    
    P = Gtr.*conj(Gt); %P coeff. eq(14)
    
    Gt2 = abs(P.^(1/3)).*exp(1i.*angle(P)); %next guess. eq(15)
    
    Gw2 = fft(Gt2,[],2).*exp(-1i.*phase); %to frequency domain and remove phase from glass. eq(16) and (17)
    
    GGuess = sum(Gw2.*intzz,1);
    
    GGuess = sqrt(Int).*exp(1i.*angle(GGuess)); %multiply by fundamental
    
    %error estimation
    
    retr = SHGv21(phase,abs(GGuess).^2,angle(GGuess)); %normalized retrieved scan
        
    Err=1-sum(sum(sqrt(img.*retr)))/sqrt(sum(sum(retr))*sum(sum(img)));
     
    iter = iter+1;
    
    figure(2);
    colormap(parula)
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
    if iter>=N
        break
    end
    
end
%% Plot
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
xlim([1.4 3.1])
legend('Retrieved', 'Measured','Phase')
hold off

%temporal domain
Np = length(wf);
% fs = wf(end)*2/2/pi;

% dt = 1/fs;

% t =(-Np/2:Np/2-1)*dt;
%     Np = length(w_exp);
% dw = abs(mean(diff(wf)));
% dt = 2*pi/dw;
% t =(-Np/2:Np/2-1)*dt;
% t = linspace(-(Np-1)*dt/Np,(Np-1)*dt/Np,Np);

dt = mean(2*pi*(wf(1)- wf(end))./wf.^2);
t = linspace(-dt*(Np-1),dt*(Np-1),Np);

figure(4);
Ret = ifft(GGuess);
Ret = ifftshift(Ret);
phase = unwrap(angle(Ret));
t2 = linspace(t(1),t(end),1920);
Inten = abs(Ret).^2./max(abs(Ret).^2);


Inten = interp1(t,Inten,t2,'pchip');
phase = interp1(t,phase,t2,'pchip');
[val,idx] = max(Inten);
t2 = t2 + t2(end-idx);
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





