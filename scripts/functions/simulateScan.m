function [scan,Ew,Et] = simulateScan(Shape,tau,phase,w0,wf,z)
%SIMULATESCAN Simulate a dscan
%Inputs:
%Shape: Pulse profile:
%                   - 'gauss' - Gaussian
%                   - 'sech' - sech^2
%tau: pulse duration in fs
%phase: phase function. can be generated e.g. by Jan's script generate_phas
%w0 - central frequency, rad/fs
%wf - frequency vector, rad/fs
%z - glass insertion range, mm
%   Detailed explanation goes here


c = 3e8;


%generate input fields
switch Shape
    case 'gauss'
        sigma = 8*log(2)/(tau^2);
%         tg  = tau/(2*log(2));
%         Ew = exp(-((wf-w0).*tg./2).^2).*exp(1i.*phase);
        Ew = exp(-((wf-w0).^2)./sigma).*exp(1i.*phase);
    case 'sech'
        ts = tau/1.763;
        Ew = sech(pi.*(wf-w0).*ts./2).^2.*exp(1i.*phase);
        
end

%phase matrix
wf = linspace(wf(1)-0.5,wf(end)+0.5,length(wf));
wl_fun = (2*pi*c./wf).*1e-9; %convert freq. to wl to calculate refr.index 
% wl_fun = linspace(wl_fun(1),wl_fun(end),length(wl_fun));
n = sqrt(1+1.03961212./(1-0.00600069867./wl_fun.^2)+0.231792344./(1-0.0200179144./wl_fun.^2)+1.01046945./(1-103.560653./wl_fun.^2)); %BK7 refractive index

k = wf.*1e15.*n./c; %wavenumber

zk = (z*1e-3)'*k; %kz matrix

pm = exp(1i.*zk);

%generate scan
Ewf = kron(Ew,ones(length(z),1)).*pm; %extend pulse vector into matrix

dw = abs(mean(diff(wf)));
Np = length(wf);
t = (dw./2*pi)*(0:length(Np)-1)/length(Np);


Et = ifft(Ewf,[],2);

Eshg = Et.^2;

dscan = fft(Eshg,[],2);
dscan = fftshift(dscan,2);
scan = abs(dscan).^2./max(max(abs(dscan).^2.));

figure();
imagesc(wf+wf(1),z,scan)

Et = ifft(Ew);
Et = ifftshift(Et);

end



