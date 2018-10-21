function [points,features]=peak_detector(t, X, varargin)
%PEAK_DETECTOR A peak detector for oscillating timeseries. The timeseries is divided into periods measured between
%successive minima. Amplitude and time window sizes can be set to select scale of peaks.
%
% If results are plotted, use left/right arrows to switch between columns
% of X.
%
% Required Inputs:
%  t - vector of timepoints
%  X - vector matching t, or matrix with columns matching t
%
% Optional input:
%  delta - fraction of amplitude for threshold: 
%     - x is max if x>last min + delta
%     scalar {default, f=0.5}
%  wsz - window size within which to check for extrema
%
% Optional name-value pair arguments:
%  'Interpolate', {true}/false - for sub-resolution peak locations
%  'ThresholdPercentiles', default=[0,100] (threshold is fraction of this range of amplitudes)
%  'Plot',true/{false} (use to control plotting with output arguments)
%
% Outputs:
%  features - struct containg feature distributions for each column of X
%  points - struct containing special points, for plotting

%check inputs and parse optional parameters
% [thrPtiles,delta,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

if isvector(X)
    X=X(:);
end

[delta,platThresh,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

nX=size(X,2); %number of traces
dt=t(2)-t(1); %TODO: plateau features assume DT is constant!

thrPtiles=[0,100];
globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

deltaIsFraction=true;
if deltaIsFraction
    delta=delta*globalXamp;
end

% hold off; plot(t,X); hold on

points(nX)=struct('tPer',[],'xPer',[],'iPer',[],...
    'tUp',[],'xUp',[],'iUp',[],'tDown',[],'xDown',[],'iDown',[],...
    'tMax',[],'xMax',[],'iMax',[],'tMin',[],'xMin',[],'iMin',[]);

% maxIsNext=zeros(1,nX); %look for a minimum first
maxIsNext=ones(1,nX); %look for a maximum first

lastMax=X(1,:); %initialize first point
lastMin=X(1,:);
maxix=ones(1,nX);
minix=ones(1,nX);
for i=2:length(t)
    this=X(i,:);
    [lastMax,maxrow]=max([this;lastMax],[],1); %update most recent maximum for all traces
    maxix(maxrow==1)=i; %maxrow(j)==1 means this(j) is new max, maxrow(j)==2 means lastMax(j) was bigger than this(j).
    [lastMin,minrow]=min([this;lastMin],[],1);
    minix(minrow==1)=i;
    
    maxFound=maxIsNext & this<lastMax-delta;
    idx=find(maxFound); %for each trace with max found, store the point
    for j=idx
        points(j).iMax(end+1)=maxix(j);
        points(j).tMax(end+1)=t(maxix(j));
        points(j).xMax(end+1)=lastMax(j);
        lastMin(maxFound)=this(maxFound);
        minix(maxFound)=i;
        maxIsNext(maxFound)=false;
%         plot(points(j).tMax(end),points(j).xMax(end),'v');
    end
    
    minFound=~maxIsNext & this>lastMin+delta;
    idx=find(minFound);
    for j=idx
        points(j).iMin(end+1)=minix(j);
        points(j).tMin(end+1)=t(minix(j));
        points(j).xMin(end+1)=lastMin(j);
        lastMax(minFound)=this(minFound);
        maxix(minFound)=i;
        maxIsNext(minFound)=true;
%         plot(points(j).tMin(end),points(j).xMin(end),'^');
    end
    
%     plot(t(maxix),lastMax,'v');
%     plot(t(minix),lastMin,'^');
end

features=struct('period',[],'APD',[],'PF',[],'amp',[],'baseline',[],'peaks',[],'pthresh',[]);
for i=1:nX
    
    if numel(points(i).tMin)>1
        
    %trim leading and trailing maxima
    if points(i).tMax(1)<points(i).tMin(1)
        points(i).iMax=points(i).iMax(2:end);
        points(i).tMax=points(i).tMax(2:end);
        points(i).xMax=points(i).xMax(2:end);
    end
    if points(i).tMax(end)>points(i).tMin(end)
        points(i).iMax=points(i).iMax(1:end-1);
        points(i).tMax=points(i).tMax(1:end-1);
        points(i).xMax=points(i).xMax(1:end-1);
    end
    
    points(i).iPer=points(i).iMin;
    points(i).tPer=points(i).tMin;
    points(i).xPer=points(i).xMin;
    
    features(i).range=globalXamp(i);
    features(i).period=diff(points(i).tMin); %period defined from minimum to minimum
    
    %simple method: baseline=average of successive minima
    features(i).baseline=0.5*(points(i).xMin(1:end-1)+points(i).xMin(2:end));
    features(i).peaks=points(i).xMax;
    features(i).amp=features(i).peaks-features(i).baseline;
    
    %interpolate active phase threshold crossing per period
    %interp of times is needed for good visual
    %local linear detrend using two minima?
    nPer=length(features(i).period);
%     features(i).active=false(size(t)); %indicators for active/silent (for plotting)
    for j=1:nPer
        ix=points(i).iMin(j):points(i).iMin(j+1);
        tt=t(ix);
        xx=X(ix,i);
        
        thr=features(i).baseline(j)+platThresh*features(i).amp(j);
        up=xx>=thr;
        iup=find(up,1,'first'); %assume any intermediate crossings still part of the up state
        idwn=find(up,1,'last');
        
        %interpolate t - note if baseline is super uneven, can be first or
        %last element: gives error, need to index into full t/X
        
%         tup=interp1(xx(iup-1:iup),tt(iup-1:iup),thr);
%         tdwn=interp1(xx(idwn:idwn+1),tt(idwn:idwn+1),thr);

%         tup=interp1(X(ix(iup-1:iup)),t(ix(iup-1:iup)),thr);
%         tdwn=interp1(X(ix(idwn:idwn+1)),t(ix(idwn:idwn+1)),thr);
        
        tup=(thr-X(ix(iup)-1,i))*(t(ix(iup))-t(ix(iup)-1))/(X(ix(iup),i)-X(ix(iup)-1,i))+t(ix(iup)-1);
        tdwn=(thr-X(ix(idwn),i))/(X(ix(idwn)+1,i)-X(ix(idwn),i))*(t(ix(idwn)+1)-t(ix(idwn)))+t(ix(idwn));

%         tup=(thr-xx(iup-1))*(tt(iup)-tt(iup-1))/(xx(iup)-xx(iup-1))+tt(iup-1);
%         tdwn=(thr-xx(idwn))*(tt(idwn+1)-tt(idwn))/(xx(idwn+1)-xx(idwn))+tt(idwn);
                
        points(i).iUp(j)=iup;
        points(i).tUp(j)=tup;
        points(i).xUp(j)=thr;
        points(i).iDown(j)=idwn;
        points(i).tDown(j)=tdwn;
        points(i).xDown(j)=thr;
        
        features(i).pthresh(j)=thr;
%         features(i).active(ix(up))=1; %needs interpolation... handle multiple crossings
    end
    
    features(i).APD=points(i).tDown-points(i).tUp;
    features(i).PF=features(i).APD./features(i).period;
   
    else
        %had less than two minima: can't compute features.
        features(i).period=0;
        features(i).baseline=0;
        features(i).peaks=0;
        features(i).amp=0;
        features(i).APD=0;
        features(i).PF=0;
    end
end



%plot to show performance
if nargout==0||doPlot==1
    tix=1;
    if ~exist('figID','var')||isempty(figID)
        figID=gcf; %new figure if none available, otherwise current fig
    else
        figID=figure(figID);
    end
    if dokeypress
        figID.KeyPressFcn=@keypressFcn;
    end
    
    plotData()
end

%nested functions can see variables in caller's scope
    function plotData()
%         clf
        delete(findobj(gca,'Tag','plateau_detector'));

        thisX=X(:,tix);
        
        plot(t,thisX,'k-','Tag','plateau_detector')
        axis tight
        hold on
        
        if length(points(tix).tMin)>1

        tupdwn=[points(tix).tUp; points(tix).tDown];
        ythresh=[1;1]*features(tix).pthresh;
        plot(tupdwn,ythresh,'r-','Tag','plateau_detector')
        
        tper=[points(tix).tMin(1:end-1);points(tix).tMin(2:end)];
        ybase=[1;1]*features(tix).baseline;
        plot(tper,ybase,'b-','Tag','plateau_detector')
        
%         %add the up/down times into a supplemented t and x for plotting
%         % note: still misses intermediate up/down crossings
%         tsup=[t',points(tix).tUp,points(tix).tDown];
%         xsup=[thisX',points(tix).xUp,points(tix).xDown];
%         [tsup,ixs]=sort(tsup);
%         xsup=xsup(ixs);
%         xAct=nan(size(tsup));
%         xSil=nan(size(tsup)); 
%         for int=1:length(features(tix).period)
%             isx=find(tsup>=points(tix).tMin(int)&tsup<=points(tix).tMin(int+1));
%             xsx=xsup(isx);
% 
%             sup=xsx>=features(tix).pthresh(int);
%             sdwn=xsx<=features(tix).pthresh(int);
%             xAct(isx(sup))=xsup(isx(sup));
%             xSil(isx(sdwn))=xsup(isx(sdwn));
%         end
%         plot(tsup,xAct,'r-','Tag','plateau_detector')
%         plot(tsup,xSil,'b-','Tag','plateau_detector')
        
        plot(points(tix).tMax,points(tix).xMax,'rv','Tag','plateau_detector')
        plot(points(tix).tMin,points(tix).xMin,'r^','Tag','plateau_detector')
%         plot(points(tix).tDXMax,points(tix).xDXMax,'g>','Tag','plateau_detector')
%         plot(points(tix).tDXMin,points(tix).xDXMin,'g<','Tag','plateau_detector')

        end
        xlabel('t')
        ylabel('x')
        YLIM=ylim();
        ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
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

function [delta,platThresh,doPlot,figID,dokeypress]=parseArgs(t, X, varargin)

    %default parameters
    defaultDelta=0.35;
    defaultPlatThresh=0.5;
    
    doPlot=false;
    doKeypress=true;
    figID=[];
    
    p=inputParser;
    validX=@(x) isreal(x) && size(x,1)==length(t); %traces are columns of X
    validSwitch=@(x) isscalar(x) && (isnumeric(x)||islogical(x));
    addRequired(p,'t',@(x) isreal(x));
    addRequired(p,'X',validX);
    addOptional(p,'delta',defaultDelta);
    addOptional(p,'platThresh',defaultPlatThresh);
    addParameter(p,'Plot',doPlot,validSwitch);
    addParameter(p,'Keypress',doKeypress,validSwitch);
    addParameter(p,'FigureID',figID);
    
    parse(p,t,X,varargin{:});
    
    delta=p.Results.delta;
    platThresh=p.Results.platThresh;
    doPlot=p.Results.Plot;
    dokeypress=p.Results.Keypress;
    figID=p.Results.FigureID;
    
end