function [retr,field] = retrtry(img,wf,z,N,fund)
%RETRFUN Summary of this function goes here
% retrieval function. based on fast retrieval algorithm by miguel.

% inputs:
% scan - processed image of measured dscan already calibrated on frequency
% axis. 
% wf - fundamental frequency vector
% z - glass insertion range
% N - maximum number of iterations in the algorithm

%outputs:
% retr  - retrieved dscan
% field - retrieved field



%% Initialization & constants

% constants

inten  = abs(fund.^2)./max(abs(fund.^2));

w3 = wf+wf(1); %shg frequencies

tau = 15; %FWHM of the first guess in time domain in fs

c = 3e8; %speed of light in vacuum, m/s

w0 = 2.35; %central frequency, rad/fs

sigma = 8*log(2)/(tau^2); %std.dev in gaussian pulse 1/fs^2

GVD = 50.621; %BK7 group velocity dispersion, fs^2/mm

wl_fun = (2*pi*c./wf).*1e-9; %convert freq. to wl to calculate refr.index 

% wl_fun = linspace(wl_fun(1),wl_fun(end),length(wl_fun));

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

int_z=sum(img,2);

int_z=int_z/sum(int_z);

intzz=wextend('ac','sp0',int_z,length(wf)-1,'d');

idz=int_z==max(int_z);

amp_trace = sqrt(img); %measured d scan field

% pm = exp(1i.*zk); %phase matrix. eq(7)

% Gaussian pulse 
GGuess = sqrt(exp(-((wf-w0).^2)./sigma)).*exp(1i.*phase); % eq(6)

% GGuess=abs(ifft(sqrt(fft(amp_trace(idz,:),[],2)),[],2));
% 
%% retrieval
%extend gauss pulse vector into matrix





while Err > 5e-5
    
    if iter >= N
        break
    end
    
    amp_i = kron(GGuess,ones(length(z),1)); % extend fund vector into a matrix
    
    U_i=ifft(amp_i.*exp(1i*zk),[],2); % multiply by phase matrix and ifft to time domain
    
    S_shg=fft(U_i.^2,[],2); %shg and go to freq domaing

    S_f=amp_trace.*exp(1i*angle(S_shg)); %multiply by scan amplitude and keep the phase
    
    U_t=ifft(S_f,[],2); %to time domain
    
    P=U_t.*conj(U_i); 
    
    UU_t=abs(P.^(1/3)).*exp(1i*angle(P));
    
    UU_f=fft(UU_t,[],2).*exp(-1i*zk);
    
    E_f2=sum(UU_f.*intzz,1);
    
    GGuess=sqrt(inten).*exp(1i*angle(E_f2));
    %error estimation
    
    retr = abs(S_shg).^2./max(max(abs(S_shg).^2)); %normalized retrieved scan
    
    img1 = img./max(max(img)); %normalized measured
    
    mu = sum(sum(img1.*retr))./sum(sum(retr))+eps; %minimization vector/constant
    
    Err = sum(sum((img1 - mu.*retr).^2))./(length(z)*length(wf)); %error estimation
    
    iter = iter+1;
    
    %plot the trace
    figure(2);
    
    colormap(parula)
    
    subplot(2,1,1)
    imagesc(w3,z,retr)
    subplot(2,1,2)
    imagesc(w3,z,img)
    title(['iter=',num2str(iter), 'Error=',num2str(Err)])
    ylabel('glass,mm')
    xlabel('frequency, rad/fs')
    drawnow;
    
    
end

%% output

field = GGuess;

    
end

