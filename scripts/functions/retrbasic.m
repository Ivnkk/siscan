function [retr,field] = retrbasic(img,wf,z,N,fund)
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




%% Initialization & constants

% constants
iter = 1;

Err = 1;

fshg = wf+min(wf); %shg frequencies

wl=300./wf; %wavelength vector. spacing?

phase=2*pi*wextend('ac','sp0',z',length(wf)-1,'r').*nBK7(wextend('ar','sp0',wl,length(z)-1,'d')/1000)./(wextend('ar','sp0',wl,length(z)-1,'d')/1000000); % kz matrix 

Meas=sqrt(img); %amplitude of the trace

int_z=sum(img,2);

int_z=int_z/sum(int_z);

intzz=wextend('ac','sp0',int_z,length(wf)-1,'d');

idz=int_z==max(int_z);

GGuess=abs(ifft(sqrt(fft(Meas(idz,:),[],2)),[],2));
%% retrieval

while Err>1e-4
    
    GGuess = kron(GGuess,ones(length(z),1)); %extend pulse vector into matrix
    
    Gt = ifft(GGuess.*exp(1i.*phase),[],2); %to time domain. eq(8) and eq(9)
    
    Gshgw = fft(Gt.^2,[],2); %generate shg and go to frequency domain. eq(11)
    
    Gup = Meas.*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)
    
    Gtr = ifft(Gup,[],2); %to time domain. eq(13)
    
    P = Gtr.*conj(Gt); %P coeff. eq(14)
    
    Gt2 = abs(P.^(1/3)).*exp(1i.*angle(P)); %next guess. eq(15)
    
    Gw2 = fft(Gt2,[],2).*exp(-1i.*phase); %to frequency domain and remove phase from glass. eq(16) and (17)
    
    GGuess = sum(Gw2.*intzz,1);
    
    GGuess = sqrt(fund).*exp(1i.*angle(GGuess));
    
    %error estimation
    
    retr = SHGv21(phase,abs(GGuess).^2,angle(GGuess)); %normalized retrieved scan
        
    Err=1-sum(sum(sqrt(img.*retr)))/sqrt(sum(sum(retr))*sum(sum(img)));
     
    iter = iter+1;
    
    figure(2);
    colormap(parula)
    subplot(2,1,1)
    imagesc(fshg,z,retr)
    subplot(2,1,2)
    imagesc(fshg,z,img)
    title(['iter=',num2str(iter), 'Error=',num2str(Err)])
    ylabel('glass,mm')
    xlabel('frequency, rad/fs')
    drawnow;
    if iter>=N
        break
    end
    
end

%% output

%     plot the results
%     figure(2);
%     colormap(parula)
%     subplot(2,1,1)
%     imagesc(fshg,z,retr)
%     subplot(2,1,2)
%     imagesc(fshg,z,img)
%     title(['iter=',num2str(iter), 'Error=',num2str(Err)])
%     ylabel('glass,mm')
%     xlabel('frequency, rad/fs')

field = GGuess;

    
end

