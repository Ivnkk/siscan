
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

addpath('G:\Atto\Data\LASC\MHz\mhzharmonics\2019-05-22\'); 
addpath('D:\PhD\Programmes\MATLAB\siscan\scripts\functions');
addpath('D:\PhD\Programmes\MATLAB\siscan\scripts\data');

%% Load the data 

load('wavelength.mat') % imaging spectrometer wavelength calibration file

load('fund3.mat') % fundamental spectrum

Int = fnspec(:,2)'; %fundamental spectrum intensity

f = fnspec(:,1)'; %fundamental spectrum frequency 

img = imread('1st1','png'); %scan to load





%% Initialization & constants

N = 20; %max number of iterations

iter = 1; % counter

Err = 1; %rms error starting value

fshg = f+min(f); %shg frequencies

wl_f = 300./f; %wavelength vector. spacing?

z = linspace(-2,2,256); % glass insertion range.

img = processScan(img,fnspec,z,wl);

img = max(0,img);

phase = 2*pi*wextend('ac','sp0',z',length(f)-1,'r').*nBK7(wextend('ar','sp0',wl_f,length(z)-1,'d')/1000)./(wextend('ar','sp0',wl_f,length(z)-1,'d')/1000000); % kz matrix 

Meas=sqrt(img); %amplitude of the trace

int_z=sum(img,2); 

int_z=int_z/sum(int_z); %normalize

intzz=wextend('ac','sp0',int_z,length(f)-1,'d');

idz=int_z==max(int_z);

GGuess=abs(ifft(sqrt(fft(Meas(idz,:),[],2)),[],2));
%% retrieval

while Err>0.001
    
    GGuess = kron(GGuess,ones(length(z),1)); %extend pulse vector into matrix
    
    Gt = ifft(GGuess.*exp(1i.*phase),[],2); %to time domain. eq(8) and eq(9)
    
    Gshgw = fft(Gt.^2,[],2); %generate shg and go to frequency domain. eq(11)
    
    Gup = Meas.*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)
    
    Gtr = ifft(Gup,[],2); %to time domain. eq(13)
    
    P = Gtr.*conj(Gt); %P coeff. eq(14)
    
    Gt2 = abs(P.^(1/3)).*exp(1i.*angle(P)); %next guess. eq(15)
    
    Gw2 = fft(Gt2,[],2).*exp(-1i.*phase); %to frequency domain and remove phase from glass. eq(16) and (17)
    
    GGuess = sum(Gw2.*intzz,1);
    
%     GGuess = sqrt(Int).*exp(1i.*angle(GGuess)); %multiply by fundamental
    
    %error estimation
    
    retr = SHGv21(phase,abs(GGuess).^2,angle(GGuess)); %normalized retrieved scan
        
    Err=1-sum(sum(sqrt(img.*retr)))/sqrt(sum(sum(retr))*sum(sum(img)));
     
    iter = iter+1;
    
    figure(2);
    colormap(parula)
    subplot(2,1,1)
    imagesc(fshg,z,retr)
    xlim([0.6 1])
    subplot(2,1,2)
    imagesc(fshg,z,img)
    title(['iter=',num2str(iter), 'Error=',num2str(Err)])
    ylabel('glass,mm')
    xlabel('frequency, PHz')
    xlim([0.6 1])
    drawnow;
    if iter>=N
        break
    end
    
end
%% Plot
phase = angle(GGuess);
phase = unwrap(phase);

figure(3);
plot(f,abs(GGuess).^2./max(abs(GGuess).^2.))
title('Spectral domain')
hold on
plot(f,Int,'k--')
yyaxis left
ylabel('Spectral power, a.u.')
xlabel('freq, PHz')
legend('Retrieved', 'Measured')
yyaxis right

plot(f,phase)
ylabel('Spectral phase, rad')
xlim([0.2 0.6])
legend('Retrieved', 'Measured','Phase')
hold off

%temporal domain
Np = length(f);
% fs = wf(end)*2/2/pi;

% dt = 1/fs;

% t =(-Np/2:Np/2-1)*dt;
%     Np = length(w_exp);
% dw = abs(mean(diff(wf)));
% dt = 2*pi/dw;
% t =(-Np/2:Np/2-1)*dt;
% t = linspace(-(Np-1)*dt/Np,(Np-1)*dt/Np,Np);

dt = 1./abs((f(1)- f(end)));
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





