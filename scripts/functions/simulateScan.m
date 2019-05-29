function [scan,Ew] = simulateScan(Shape,tau,phase,w0,wf,z)
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



f0 = w0/2/pi;
f = wf/2/pi;

%generate input fields
switch Shape
    case 'gauss'
%         sigma = 8*log(2)/(tau^2);
        tg  = tau/(2*log(2));
        Ew = sqrt(exp(-(2*pi*(f-f0).*tg/2).^2)).*exp(1i.*phase);
%         Ew = sqrt(exp(-((wf-w0).^2)./sigma)).*exp(1i.*phase);
    case 'sech'
        ts = tau/1.763;
        Ew = sqrt(sech(pi.*(f-f0).*ts./2).^2).*exp(1i.*phase);
        
end

%phase matrix
f=linspace(min(f)-min(f)/2,max(f) + max(f/2),length(f)); % extend the frequency vector to make sure that the 
wl_fun = 300./f; %convert freq. to wl to calculate refr.index 

zk=2*pi*wextend('ac','sp0',z',length(f)-1,'r').*nBK7(wextend('ar','sp0',wl_fun,length(z)-1,'d')/1000)./(wextend('ar','sp0',wl_fun,length(z)-1,'d')/1000000); % phase matrix
pm = exp(1i.*zk);

%generate scan
Ewf = kron(Ew,ones(length(z),1)).*pm; %extend pulse vector into matrix

Et = ifft(Ewf,[],2);

Eshg = Et.^2;

dscan = fft(Eshg,[],2);
dscan = fftshift(dscan,2);
scan = abs(dscan).^2./max(max(abs(dscan).^2.));

figure();
imagesc(f+min(f),z,scan)

end



