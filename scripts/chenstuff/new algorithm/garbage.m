

% tmax = 40;
% t = linspace(-5.*tmax,5.*tmax,1024);
% N  = length(t);% number of points
% dt = mean(diff(t));% dt
% 
% df = 0.5./tmax;
% 
% w = linspace(0,df*N,N);
% 
% x = sin(2*pi*w(5).*t);
% 
% xw = fft(x);
% plot(w,abs(xw))
% 

t = 0:1/50:10-1/50;                     
x = sin(2*pi*15*t) + sin(2*pi*20*t);
plot(t,x)