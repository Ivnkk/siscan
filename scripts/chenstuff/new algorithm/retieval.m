clear;
clc;
close all;

load('sample0.mat'); %load the simulated scan

wl=300./f0; %wavelength vector. spacing?

dphase=2*pi*wextend('ac','sp0',glass',length(f0)-1,'r').*nBK7(wextend('ar','sp0',wl,length(glass)-1,'d')/1000)./(wextend('ar','sp0',wl,length(glass)-1,'d')/1000000); % kz matrix 

amp_trace=sqrt(intshg); %amplitude of the trace

int_z=sum(intshg,2);

int_z=int_z/sum(int_z);

intzz=wextend('ac','sp0',int_z,length(f0)-1,'d');

idz=int_z==max(int_z);

E_f=abs(ifft(sqrt(fft(amp_trace(idz,:),[],2)),[],2));

figure(1)
subplot(2,1,1)
imagesc(fshg,glass,intshg);

for i=1:1000

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
       
    shg_r=SHGv21(dphase,abs(E_f).^2,angle(E_f));
    
    error=1-sum(sum(sqrt(intshg.*shg_r)))/sqrt(sum(sum(shg_r))*sum(sum(intshg)));
    
    subplot(2,1,2)
    
    imagesc(fshg,glass,shg_r);
    
    title(['i=',num2str(i),' error=',num2str(error)]);

    drawnow;
    
    if error<5e-5
        break;
    end
end
    