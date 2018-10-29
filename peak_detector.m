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

[delta,platThresh,minAmp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

nX=size(X,2); %number of traces
dt=t(2)-t(1); %TODO: plateau features assume DT is constant!

thrPtiles=[0,100];
globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

sufficientGlobalAmp=globalXamp>minAmp;

deltaIsFraction=true;
if deltaIsFraction
    delta=delta*globalXamp;
end

% hold off; plot(t,X); hold on
DX=slopeY(t,X); %todo: noise-robust method; measure between tmin(i) and tmin(i+1)?

pt=struct('ix',[],'t',[],'x',[],'dx',[]);

points=repmat(struct('period',pt,'up',pt,'down',pt,...
    'max',pt,'min',pt,'dxmax',pt,'dxmin',pt),1,nX);
% maxIsNext=zeros(1,nX); %look for a minimum first (often calls first point
maxIsNext=ones(1,nX); %look for a maximum first (most reliable, but appears to miss some cases)

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
    
    maxFound=maxIsNext & this<lastMax-delta & sufficientGlobalAmp;
    minFound=~maxIsNext & this>lastMin+delta & sufficientGlobalAmp;
    
    idx=find(maxFound); %for each trace with max found, store the point
    for j=idx
        points(j).max.ix(end+1)=maxix(j);
        points(j).max.t(end+1)=t(maxix(j));
        points(j).max.x(end+1)=lastMax(j);
        lastMin(maxFound)=this(maxFound);
        minix(maxFound)=i;
        maxIsNext(maxFound)=false;
%         plot(points(j).max.t(end),points(j).max.x(end),'v');
    end
    
    idx=find(minFound);
    for j=idx
        points(j).min.ix(end+1)=minix(j);
        points(j).min.t(end+1)=t(minix(j));
        points(j).min.x(end+1)=lastMin(j);
        lastMax(minFound)=this(minFound);
        maxix(minFound)=i;
        maxIsNext(minFound)=true;
%         plot(points(j).min.t(end),points(j).min.x(end),'^');
    end
    
%     plot(t(maxix),lastMax,'v');
%     plot(t(minix),lastMin,'^');
end


features=struct('period',[],'APD',[],'PF',[],'amp',[],'baseline',[],'peaks',[],'pthresh',[],'maxslope',[],'minslope',[]);
for i=1:nX
        
    %trim leading and trailing maxima
    if  points(i).max.t(1)<points(i).min.t(1)
        points(i).max.ix=points(i).max.ix(2:end);
        points(i).max.t=points(i).max.t(2:end);
        points(i).max.x=points(i).max.x(2:end);
    end
    if points(i).max.t(end)>points(i).min.t(end)
        points(i).max.ix=points(i).max.ix(1:end-1);
        points(i).max.t=points(i).max.t(1:end-1);
        points(i).max.x=points(i).max.x(1:end-1);
    end
    
    points(i).period.ix=points(i).min.ix;
    points(i).period.t=points(i).min.t;
    points(i).period.x=points(i).min.x;
    nPer=length(points(i).period.t)-1;

    if nPer>=1
    
    features(i).period=diff(points(i).period.t); %period defined from minimum to minimum
    
    %simple method: baseline=average of successive minima
    features(i).baseline=0.5*(points(i).min.x(1:end-1)+points(i).min.x(2:end));
    features(i).peaks=points(i).max.x;
    features(i).amp=features(i).peaks-features(i).baseline;
    
    %interpolate active phase threshold crossing per period
    %interp of times is needed for good visual
    %local linear detrend using two minima?
%     features(i).active=false(size(t)); %indicators for active/silent (for plotting)
    for j=1:nPer
        ix=points(i).min.ix(j):points(i).min.ix(j+1);
        tt=t(ix);
        xx=X(ix,i);
        dx=DX(ix,i);
        
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
        
        %direct linear interpolation
        tup=(thr-X(ix(iup)-1,i))*(t(ix(iup))-t(ix(iup)-1))/(X(ix(iup),i)-X(ix(iup)-1,i))+t(ix(iup)-1);
        tdwn=(thr-X(ix(idwn),i))/(X(ix(idwn)+1,i)-X(ix(idwn),i))*(t(ix(idwn)+1)-t(ix(idwn)))+t(ix(idwn));

%         tup=(thr-xx(iup-1))*(tt(iup)-tt(iup-1))/(xx(iup)-xx(iup-1))+tt(iup-1);
%         tdwn=(thr-xx(idwn))*(tt(idwn+1)-tt(idwn))/(xx(idwn+1)-xx(idwn))+tt(idwn);
                
        points(i).up.ix(j)=iup;
        points(i).up.t(j)=tup;
        points(i).up.x(j)=thr;
        points(i).down.ix(j)=idwn;
        points(i).down.t(j)=tdwn;
        points(i).down.x(j)=thr;
        
        features(i).pthresh(j)=thr;
%         features(i).active(ix(up))=1; %needs interpolation... handle multiple crossings


        [dxmax,idxmax]=max(dx);
        points(i).dxmax.ix(j)=ix(1)+idxmax-1;
        points(i).dxmax.t(j)=t(ix(idxmax));
        points(i).dxmax.x(j)=X(ix(idxmax),i);
        points(i).dxmax.dx(j)=dxmax;
        
        [dxmin,idxmin]=min(dx);
        points(i).dxmin.ix(j)=ix(1)+idxmin-1;
        points(i).dxmin.t(j)=t(ix(idxmin));
        points(i).dxmin.x(j)=X(ix(idxmin),i);
        points(i).dxmin.dx(j)=dxmin;
    end
    
    features(i).APD=points(i).down.t-points(i).up.t;
    features(i).PF=features(i).APD./features(i).period;
   
    features(i).range=globalXamp(i);
    features(i).maxslope=points(i).dxmax.dx;
    features(i).minslope=points(i).dxmin.dx;
    else
        %had less than two minima: can't compute features.
        features(i).period=[];
        features(i).baseline=[];
        features(i).peaks=[];
        features(i).amp=[];
        features(i).APD=[];
        features(i).PF=[];
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
        
        if length(points(tix).min.t)>1

        tupdwn=[points(tix).up.t; points(tix).down.t];
        ythresh=[1;1]*features(tix).pthresh;
        plot(tupdwn,ythresh,'r-','Tag','plateau_detector')
        
%         tper=[points(tix).min.t(1:end-1); points(tix).min.t(2:end)];
%         ybase=[1;1]*features(tix).baseline;
%         plot(tper,ybase,'b-','Tag','plateau_detector')
        
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
        
        plot(points(tix).max.t,points(tix).max.x,'rv','Tag','plateau_detector')
        plot(points(tix).min.t,points(tix).min.x,'r^','Tag','plateau_detector')
        plot(points(tix).dxmax.t,points(tix).dxmax.x,'g>','Tag','plateau_detector')
        plot(points(tix).dxmin.t,points(tix).dxmin.x,'g<','Tag','plateau_detector')

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

function [delta,platThresh,minAmp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin)

    %default parameters
    defaultDelta=0.35;
    defaultPlatThresh=0.5;
    defaultminAmp=0;
    
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
    addParameter(p,'MinimumAmplitude',defaultminAmp,@isreal);
    addParameter(p,'Plot',doPlot,validSwitch);
    addParameter(p,'Keypress',doKeypress,validSwitch);
    addParameter(p,'FigureID',figID);
    
    parse(p,t,X,varargin{:});
    
    delta=p.Results.delta;
    platThresh=p.Results.platThresh;
    minAmp=p.Results.MinimumAmplitude;
    doPlot=p.Results.Plot;
    dokeypress=p.Results.Keypress;
    figID=p.Results.FigureID;
    
end