function [retr,field] = retrfun(img,fund,z,N,bgvalue,use_fund)
%RETRFUN Summary of this function goes here
% retrieval function. based on fast retrieval algorithm by miguel.
% inputs:
% scan - processed 300x480 image of measured dscan already calibrated on frequency
% axis. can be created from raw 1200x1920 camera images by processScan
% function
% fund - measured fundamental spectrum - necessary for the freq. axis!
% z - glass insertion range
% bg - background value to subtract
% N - maximum number of iterations in the algorithm
% use_fund (optional) - use of the fundamental to calibrate the trace
%   Detailed explanation goes here


%% Initialization & constants
% bg subtraction

img = img-bgvalue;
img = max(0,img);



% constants

wf = fund(:,1)'; %uploaded fundamental frequency

Int = fund(:,2)'; %uploaded fundamental intensity

w3 = wf+wf(end); %shg frequencies

[~,idx] = max(Int); % find the peak of the spectrum (central frequency)

tau = 20; %FWHM of the first guess in time domain in fs

c = 3e8; %speed of light in vacuum, m/s

w0 = wf(idx); %central frequency, rad/fs

sigma = 8*log(2)/(tau^2); %std.dev in gaussian pulse 1/fs^2

GVD = 50.621; %BK7 group velocity dispersion, fs^2/mm

wl_fun = (2*pi*c./wf).*1e-9; %convert freq. to wl to calculate refr.index 

wl_fun = linspace(wl_fun(1),wl_fun(end),length(wl_fun));

n = sqrt(1+1.03961212./(1-0.00600069867./wl_fun.^2)+0.231792344./(1-0.0200179144./wl_fun.^2)+1.01046945./(1-103.560653./wl_fun.^2)); %BK7 refractive index

k = wf.*1e15.*n./c; %wavenumber

zk = (z*1e-3)'*k; %kz matrix

iter = 1; % counter

Err = 1; %initial value for error

%% First guess

%spectral phase
[~,idx] = max(img); 

GDD = GVD.*z(end-idx); 

GDD = interp1(wf,GDD,wf,'linear');

GD = cumtrapz(wf,GDD);

phase = cumtrapz(wf,GD);

% Gaussian pulse 
Gauss = exp(-((wf-w0).^2)./sigma).*exp(1i.*phase); % eq(6)

%% retrieval
GGuess = kron(Gauss,ones(length(z),1)); %extend gauss pulse vector into matrix

pm = exp(1i.*zk); %phase matrix. eq(7)

Meas = sqrt(double(img)); %measured d scan field

while Err > 0.01
    
    if iter >= N
        break
    end
    
    GGuess = GGuess.*pm; %Initial guess multiplied by pm. eq(8)
    
%     %frequency shift (if needed)
%     dw = abs(mean(diff(wf)));
%     Np = length(w_exp);
%     t = linspace(0,(Np-1)*2*pi/dw/Np,Np);
    
    Gt = ifft(GGuess,[],2); %to time domain. eq(9)
    
    Gshgt = Gt.^2; %generate shg. eq(10)
    
    Gshgw = fft(Gshgt,[],2); %to frequency domain. eq(11)
    
    Gup = Meas.*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)
    
    Gtr = ifft(Gup,[],2); %to time domain. eq(13)
    
    P = Gtr.*conj(Gt); %P coeff. eq(14)
    
    Gt2 = nthroot(abs(P),3).*exp(1i.*angle(P)); %next guess. eq(15)
    
    Gw2 = fft(Gt2,[],2).*conj(pm); %to frequency domain and remove phase from glass. eq(16) and (17)
    
    GGuess = sum(Gw2.*(z(2)-z(1)))./(z(end)-z(1)); %new field. eq(18)
    
    if use_fund == 1
        
        GGuess = sqrt(Int).*exp(1i.*angle(GGuess)); %multiply by fundamental
        
    end
    
    %error estimation
    
    retr = abs(Gshgw).^2./max(max(abs(Gshgw).^2)); %normalized retrieved scan
    
    img1 = img./max(max(img)); %normalized measured
    
    mu = sum(sum(img1.*retr))./sum(sum(retr))+eps; %minimization vector/constant
    
    Err = sum(sum((img1 - mu.*retr).^2))./(length(z)*length(wf)); %error estimation
    
    Meas = Meas./mu;
    
    iter = iter+1;
    %plot the trace
    figure(2);
    colormap(parula)
    imagesc(w3,z,retr)
    ylabel('glass,mm')
    xlabel('frequency, rad/fs')
    drawnow;
    
    
end

%% output

field = GGuess;

    
end

