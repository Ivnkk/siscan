function [retr,field] = retrbasic2(intshg,f0,glass,N,inten)
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
iter = 1;

fshg = f0+min(f0); %shg frequencies

wl=300./f0; %wavelength vector. spacing?

dphase=2*pi*wextend('ac','sp0',glass',length(f0)-1,'r').*nBK7(wextend('ar','sp0',wl,length(glass)-1,'d')/1000)./(wextend('ar','sp0',wl,length(glass)-1,'d')/1000000); % kz matrix 

amp_trace=sqrt(intshg); %amplitude of the trace

int_z=sum(intshg,2);

int_z=int_z/sum(int_z);

intzz=wextend('ac','sp0',int_z,length(f0)-1,'d');

idz=int_z==max(int_z);

E_f=abs(ifft(sqrt(fft(amp_trace(idz,:),[],2)),[],2));
%% retrieval

for i = 1:N
    
%     GGuess = kron(GGuess,ones(length(z),1)); %extend pulse vector into matrix
%     
%     Gt = ifft(GGuess.*exp(1i.*zk),[],2); %to time domain. eq(8) and eq(9)
%     
%     Gshgw = fft(Gt.^2,[],2); %generate shg and go to frequency domain. eq(11)
%     
%     Gup = Meas.*exp(1i.*angle(Gshgw)); %multiply by measured intensity. eq(12)
%     
%     Gtr = ifft(Gup,[],2); %to time domain. eq(13)
%     
%     P = Gtr.*conj(Gt); %P coeff. eq(14)
%     
%     Gt2 = abs(P.^(1/3)).*exp(1i.*angle(P)); %next guess. eq(15)
%     
%     Gw2 = fft(Gt2,[],2).*exp(-1i.*zk); %to frequency domain and remove phase from glass. eq(16) and (17)
%     
%     GGuess = sum(Gw2.*intzz,1);
% 
%     GGuess = sqrt(fund).*exp(1i.*angle(GGuess)); %multiply by fundamental
%     
%     %error estimation
%     
%     retr = SHGv21(kz,abs(GGuess).^2,angle(GGuess)); %normalized retrieved scan
%         
%     Err=1-sum(sum(sqrt(img.*retr)))/sqrt(sum(sum(retr))*sum(sum(img)));
%     
%     iter = iter+1;
%     
%     if Err < 5e-5
%         break
%     end

    amp_i = kron(E_f,ones(length(glass),1)); % extend fund vector into a matrix
    
    U_i=ifft(amp_i.*exp(1i*dphase),[],2); % multiply by phase matrix and ifft to time domain
    
    S_shg=fft(U_i.^2,[],2); %shg and go to freq domaing

    S_f=amp_trace.*exp(1i*angle(S_shg)); %multiply by scan amplitude and keep the phase
    
    U_t=ifft(S_f,[],2); %to time domain
    
    P=U_t.*conj(U_i); 
    
    UU_t=abs(P.^(1/3)).*exp(1i*angle(P));
    
    UU_f=fft(UU_t,[],2).*exp(-1i*dphase);
    
    E_f2=sum(UU_f.*intzz,1);
    
    E_f=sqrt(inten).*exp(1i*angle(E_f2));
       
    retr=SHGv21(dphase,abs(E_f).^2,angle(E_f));
    
    Err=1-sum(sum(sqrt(intshg.*retr)))/sqrt(sum(sum(retr))*sum(sum(intshg)));

    if Err<5e-5
        break;
    end
    iter = iter+1;
    
    %plot the trace
    figure(2);
    
    colormap(parula)
    
    subplot(2,1,1)
    imagesc(fshg,glass,retr)
    title(['iter=',num2str(iter), 'Error=',num2str(Err)])
    subplot(2,1,2)
    imagesc(fshg,glass,intshg)
    
    ylabel('glass,mm')
    xlabel('frequency, rad/fs')
    drawnow;
    
    
end

%% output

field = E_f;

    
end

