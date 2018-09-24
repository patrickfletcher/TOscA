function [Xdetrend,XTrend]=detrend(t,X,params)

nX=size(X,2);

dt=mode(diff(t));
fs=1/dt; 

% steepness=params.filter.steepness;
% stopAtten=params.filter.stopAtten;

XTrend=zeros(size(X));
switch params.trend.method
    case {'none'}

    case {'poly','polynomial'}
        for i=1:nX
            XTrend(:,i)=polyval(polyfit(t,X(:,i),params.trend.degree));
        end

    case {'lowpass'}
        XTrend=lowpass(X,params.trend.fpass,fs,'ImpulseResponse','iir');
%         XTrend=lowpass(X,params.trend.fpass,fs,'ImpulseResponse','iir','Steepness',steepness,'StopbandAttenuation',stopAtten);

    case {'movmean','movmedian','gaussian','lowess','loess','rlowess','rloess','sgolay'}
        wsz=round(params.trend.wsz/dt);
        XTrend=smoothdata(X,params.trend.method,wsz);
end

Xdetrend=X-XTrend;