%%% this function imports the data from the text file and reshapes the
%%% spectrum to match the dimensions of the single shot dscan
%%% INPUT: res - sensor resolution
%%%        bin - binning factor
%%%        OPTIONAL Lext/Hext - low/high extension of frequency range, IN RAD/fs!!
%%%        
%%% OUTPUT: spec:
%%%         1st column: frequency in rad/fs
%%%         2nd column: normalized intensity

%res,bin,Lext,Hext


function spec = getFundamental(res,bin,Lext,Hext) 
if nargin <3
    Lext = 0; %default values - no extension of the freq.range, just use the fundamental.
    Hext = 0; 
end
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
%% process
  c = 3e8; % speed of light
  
%   wf_1 = 4 + Hext; %border of the defined freq. range #1 

%   wf_2 = (2*pi*(c/M(end,1))*1e-6) - Lext; %border of the defined freq. range #2

  w_exp = 2*pi*(c./M(:,1))*1e-6; %measured frequency range
  
%   w = fliplr(linspace(wf_1,wf_2,length(M(:,1)))'); %axis to interpolate to
  
  w = linspace(w_exp(end), w_exp(1),length(M(:,1)))';
  
  I = fliplr(max(0,M(:,2))); % reject negative intensity values
  
  I = interp1(w,I,w,'linear')./w.^2; %reinterpolate
  
  I = I./max(I); %normalize
  
 %% resize 
 
  del = length(I) - res;
  
  w(del:end) = [];% resizing
  
  I(1:del) = [];
  
  w = arrayfun(@(i) mean(w(i:i+bin-1)),1:bin:length(w)-bin+1); %binning
  
  I = arrayfun(@(i) mean(I(i:i+bin-1)),1:bin:length(I)-bin+1);
  
  spec = [w;I]';
  
end
