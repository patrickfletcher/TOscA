function XFilt=filterTraces(t,X,method,methodparam,doPlot)

dt=mode(diff(t));
fs=1/dt; 

if ~exist('doPlot','var')
    doPlot=0;
end

switch lower(method)
    case {'none'}
        XFilt=X;
        
    case {'lowpass'}
        
        %methodparam=lowpass cutoff frequency
        if ~exist('methodparam','var')
            error('lowpass trendline requires a cutoff frequency');
        else
            fpass=methodparam;
        end
        
        XFilt=lowpass(X,fpass,fs,'ImpulseResponse','iir');
%         XFilt=lowpass(X,fpass,fs,'ImpulseResponse','iir','Steepness',steepness,'StopbandAttenuation',stopAtten);
    
    case {'movmean','movmedian','gaussian','lowess','loess','rlowess','rloess','sgolay'}
        
        %methodparam=window width
        if ~exist('methodparam','var')
            error('smoothdata trendline requires a window duration (in time units)');
        else
            wwidth=methodparam;
        end
        
        wsz=round(wwidth/dt);
        XFilt=smoothdata(X,method,wsz);
        
    case {'bandpass'}
        %TODO: edge effect corrections
        
        %methodparam=lowpass cutoff frequency
        if ~exist('methodparam','var')
            error('lowpass trendline requires a cutoff frequency');
        else
            fpass=methodparam;
        end
        
        XFilt=bandpass(X,fpass,fs,'ImpulseResponse','iir');
%         XFilt=bandpass(X,fpass,fs,'ImpulseResponse','iir','Steepness',steepness,'StopbandAttenuation',stopAtten);
end



%plot to show result
if nargout==0 || doPlot==1
    
nX=size(X,2);
tix=1;
figure('KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
%         subplot(2,1,1)
%         plot(t,X(:,tix),t,XFilt(:,tix))
%         grid on
%         xlabel('Time')
% %         ylabel('raw')
%         axis tight
%         
%         subplot(2,1,2)
%         plot(t,XFilt(:,tix))
%         grid on
%         xlabel('Time')
%         ylabel('filtered')
%         axis tight
        
        
        plot(t,X(:,tix))
        hold on
        plot(t,XFilt(:,tix),'linewidth',1.5)
        hold off
        grid on
        xlabel('Time')
        ylabel('filtered')
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