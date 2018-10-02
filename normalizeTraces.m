function XNorm=normalizeTraces(t,X,method,methodpar,doPlot)
%columns are timeseries

%
%     'zscore' - (default) normalizes by centering the data to have mean 0 
%                and scaling it to have standard deviation 1.
%
%     'center' - normalizes by centering the data to have mean 0.
%
%     'range'  - normalizes by rescaling the range of the data to the 
%                interval [0,1].

if ~exist('method','var')
    method='zscore';
    methodpar=[];
end

if ~exist('doPlot','var')
    doPlot=0;
end

switch method
    case {'none'}
        XNorm=X;
        
    case {'zscore'}
        XNorm=(X-mean(X,1))./std(X,[],1);
        
    case {'unit'}
        XNorm=(X-min(X,[],1))./(max(X,[],1)-min(X,[],1));
        if exist('methodpar','var')&&length(methodpar)==2
            intrvl=methodpar;
            XNorm=intrvl(1)+diff(intrvl)*XNorm;
        end
        
        
    case {'devmean'}
        XNorm=(X-mean(X,1))./mean(X,1);
        
    case {'devmeanpow'}
        XNorm=(X-mean(X,1))./mean(X,1).^methodpar;
        
    case {'devmedian'}
        XNorm=(X-median(X,1))./median(X,1);
        
    case {'devmedianpow'}
        XNorm=(X-median(X,1))./median(X,1).^methodpar;
        
    case {'devtrend','trend'}
        
        %no errorchecking for now..
        trendmethod=methodpar{1};
        trendpar=methodpar{2};
        
        [~,XTrend]=detrendTraces(t,X,trendmethod,trendpar);
        XNorm=(X-XTrend)./XTrend;
        
    case {'devtrend2','trend2'}
        
        %no errorchecking for now..
        trendmethod=methodpar{1};
        trendpar=methodpar{2};
        
        [~,XTrend]=detrendTraces(t,X,trendmethod,trendpar);
        XNorm=(X-XTrend)./XTrend.^2;
        
    case {'center'}
        switch methodpar
            case {'mean'}
                XNorm=X-mean(X,1);
        
            case {'median'}
                XNorm=X-median(X,1);
        
            case {'max'}
                XNorm=X-max(X,[],1);
        
            case {'min'}
                XNorm=X-min(X,[],1);
        
            otherwise
                if isnumeric(methodpar) &&isscalar(methodpar)
                    k=methodpar;
                else
                    k=0;
                end 
                XNorm=X-k;
        end
        
    case {'scale'}
        switch methodpar
            case {'mean'}
                XNorm=X./mean(X,1);
        
            case {'stdev'}
                XNorm=XNorm./std(X,[],1);
        
            case {'median'}
                XNorm=X./median(X,1);
        
            case {'max'}
                XNorm=X./max(X,[],1);
        
            case {'min'}
                XNorm=X./min(X,[],1);
        
            case {'range'}
                XNorm=XNorm./(max(X,[],1)-min(X,[],1));

            case {'iqr'}
                XNorm=XNorm./(prctile(X,25,1)-prctile(X,75,1));
        
            otherwise
                if isnumeric(methodpar) &&isscalar(methodpar)
                    k=methodpar;
                else
                    k=0;
                end 
                XNorm=X/k;
        end
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
        subplot(2,1,1)
        plot(t,X(:,tix))
        if exist('XTrend','var')
            hold on
            plot(t,XTrend(:,tix))
            hold off
        end
        grid on
        xlabel('Time')
%         ylabel('raw')
        axis tight
        
        subplot(2,1,2)
        plot(t,XNorm(:,tix))
        grid on
        xlabel('Time')
        ylabel('normalized')
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