function [time,X,opts]=loadData(filename)

if ~exist('filename','var')&&~isempty(filename)
    [filename,path]=uigetfile({'*.xls*'});
    filename=[path,filename];
end
            
[data,header]=xlsread(filename);
time=data(:,1);
X=data(:,2:end);

%decode header info
opts=[]; %placeholder

%interpolate any missing data using neighboring time points
rowIsNan=any(isnan(X),2);
timeIsNan=isnan(time);

if any(rowIsNan)
%                 flag=[flag,{'some rows have NaN'}];
    r=find(rowIsNan);
    
    if r(1)==1
        time=time(2:end,1);
        X=X(2:end,:);
        r=r(2:end)-1;
    end
    for ii=1:length(r)
        for jj=1:size(X,2)
            iix=[r(ii)-1,r(ii)+1];
            x=time(iix,1);
            v=X(iix,jj);
            xq=time(r(ii),1);
            X(r(ii),jj)=interp1(x,v,xq);
        end
    end
end