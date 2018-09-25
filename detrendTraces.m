function [Xdetrend,XTrend]=detrendTraces(t,X,method,methodparam,doPlot)

nX=size(X,2);

dt=mode(diff(t));
fs=1/dt; 

if ~exist('doPlot','var')
    doPlot=0;
end

XTrend=zeros(size(X));
switch lower(method)
    case {'none'}
%         XTrend=zeros

    case {'poly','polynomial'}
        
        %methodparam=degree of polynomial
        if ~exist('methodparam','var')||isempty(methodparam)
            degree=1;
        else
            degree=methodparam;
        end
        
        for i=1:nX
            XTrend(:,i)=polyval(polyfit(t,X(:,i),degree),t);
        end

    case {'lowpass'}
        
        %methodparam=lowpass cutoff frequency
        if ~exist('methodparam','var')||isempty(methodparam)
            error('lowpass trendline requires a cutoff frequency');
        else
            fpass=methodparam;
        end
        
        XTrend=lowpass(X,fpass,fs,'ImpulseResponse','iir');
%         XTrend=lowpass(X,params.trend.fpass,fs,'ImpulseResponse','iir','Steepness',steepness,'StopbandAttenuation',stopAtten);

    case {'movmean','movmedian','gaussian','lowess','loess','rlowess','rloess','sgolay'}
        
        %methodparam=window width
        if ~exist('methodparam','var')||isempty(methodparam)
            error('smoothdata trendline requires a window duration (in time units)');
        else
            wwidth=methodparam;
        end
        
        wsz=round(wwidth/dt);
        XTrend=smoothdata(X,method,wsz);
end

Xdetrend=X-XTrend;



%plot to show result
if nargout==0 || doPlot==1
    
nX=size(X,2);
tix=1;
figure('KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        subplot(2,1,1)
        plot(t,X(:,tix),t,XTrend(:,tix))
        grid on
        xlabel('Time')
%         ylabel('raw')
        axis tight
        
        subplot(2,1,2)
        plot(t,Xdetrend(:,tix))
        grid on
        xlabel('Time')
        ylabel('detrend')
        axis tight
        
    end

    function keypressFcn(~,event)
        switch(event.Key)
            case {'leftarrow'}
                if tix>1
                    tix=tix-1;
                    plotData()
                end
            case {'rightarrow'}
                if tix<nX
                    tix=tix+1;
                    plotData()
                end
        end
        
    end

end