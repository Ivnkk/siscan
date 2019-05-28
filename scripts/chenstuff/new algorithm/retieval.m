clear;
clc;
close all;

load('sample0.mat');
load('colorm.mat');
wl=300./f0;
dphase=2*pi*wextend('ac','sp0',glass',length(f0)-1,'r').*nBK7(wextend('ar','sp0',wl,length(glass)-1,'d')/1000)./(wextend('ar','sp0',wl,length(glass)-1,'d')/1000000);
amp_trace=sqrt(intshg);
a0=max(max(amp_trace));
amask=wextend('ar','sp0',f0*3,length(glass)-1,'d');
% amp_trace=amp_trace.*amask+0.01*a0*rand(size(amp_trace));
int_z=sum(intshg,2);
int_z=int_z/sum(int_z);
intzz=wextend('ac','sp0',int_z,length(f0)-1,'d');
idz=int_z==max(int_z);
E_f=abs(ifft(sqrt(fft(amp_trace(idz,:),[],2)),[],2));
% E_f=sqrt(inten);
% 
figure(1)
subplot(3,1,1)
imagesc(fshg,glass,intshg);
colormap(cmap0);

for i=1:1000
    amp_i=wextend('ar','sp0',E_f,length(glass)-1,'d');
    U_i=ifft(amp_i.*exp(1i*dphase),[],2);
    S_shg=fft(U_i.^2,[],2);
%     S_shg=SHGv22(dphase,E_f);
    S_f=amp_trace.*exp(1i*angle(S_shg));
    U_t=ifft(S_f,[],2);
    P=U_t.*conj(U_i);
    UU_t=abs(P.^(1/3)).*exp(1i*angle(P));
    UU_f=fft(UU_t,[],2).*exp(-1i*dphase);
    
%     E_f=mean(UU_f,1);
    E_f2=sum(UU_f.*intzz,1);
    E_f=sqrt(inten).*exp(1i*angle(E_f2));
    
    phasei=unwrap(angle(E_f));
    shg_r=SHGv21(dphase,abs(E_f).^2,angle(E_f));
    error=1-sum(sum(sqrt(intshg.*shg_r)))/sqrt(sum(sum(shg_r))*sum(sum(intshg)));
    subplot(3,1,2)
    imagesc(fshg,glass,shg_r);
    title(['i=',num2str(i),' error=',num2str(error)]);
    subplot(3,1,3);
    [AX,H1,H2]=plotyy(f0,sqrt(inten),f0,phase0-phase0(117));
    ylim(AX(2),[-10,20]);
    hold on;
    [AX,H1,H2]=plotyy(f0,abs(E_f)/max(abs(E_f)),f0,phasei-phasei(117));
    ylim(AX(2),[-10,20]);
    hold off;
    drawnow;
    
    if error<5e-5
        break;
    end
end
    