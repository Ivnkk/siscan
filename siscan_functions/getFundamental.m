%%% this function imports the data from the text file and reshapes the
%%% spectrum to match the dimensions of the single shot dscan
%%% OUTPUT: spec:
%%%         1st column: frequency in rad/fs
%%%         2nd column: normalized intensity

%res,bin,Lext,Hext


function spec = getFundamental() 

%% open the dialog and import the data
filter = {'*.dat';'*.txt'};
[file,path] = uigetfile(filter);
  if isequal(file,0)
    disp('Cancelled')
  else
    disp(['User selected ',fullfile(path,file)])
  end
  M = readtable(fullfile(path,file),'Delimiter','tab','ReadVariableNames',false);
  M = str2double(table2array(M))./1e3;
  
  wl = M(:,1)';
  wl = fliplr(wl);
  
  Int = M(:,2)';
  Int = fliplr(Int);
  Int(1774:end) = 0; %kill measured shg
  
%% process

  f = linspace(0.1,1.1,length(wl)); %desired frequency range
   
  f0 = 300./wl; % frequency range from the spectrometer
   
  I = max(0,Int); % reject negative intensity values
  
  I = interp1(f0,I,f,'linear',0)./f.^2; %reinterpolate
  
  I = I./max(I); %normalize
  
 %% resize 
 
  f(1:389) = [];% resizing
  
  f(end-399:end) = [];
  
  I(1:389) = [];% resizing
  
  I(end-399:end) = [];
 
  bin = 4;
  
  f = arrayfun(@(i) mean(f(i:i+bin-1)),1:bin:length(f)-bin+1); %binning
  
  I = arrayfun(@(i) mean(I(i:i+bin-1)),1:bin:length(I)-bin+1);
  
  spec = [f;I]';
  
end
