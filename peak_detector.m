function [F,Fdist,points,fcns]=peak_detector(t, X, varargin)
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

%TODO: input name-value support
%TODO: verify best initialization of maxIsNext toggle
%TODO: interpolate peaks/troughs?
%TODO: noise-robust method for numerical slope
%TODO: scalar feature summary output

if isvector(X)
    X=X(:);
end

[delta,platThresh,minAmp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

nX=size(X,2); %number of traces

DX=slopeY(t,X); 

thrPtiles=[0,100];
globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

sufficientGlobalAmp=globalXamp>minAmp;

deltaIsFraction=true;
if deltaIsFraction
    delta=delta*globalXamp;
end

pt=struct('ix',[],'t',[],'x',[],'dx',[]);

points=repmat(struct('period',pt,'up',pt,'down',pt,...
    'max',pt,'min',pt,'dxmax',pt,'dxmin',pt),1,nX);

%initialization of maxIsNext toggle
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

for i=1:nX
    points(i).period.ix=points(i).min.ix;
    points(i).period.t=points(i).min.t;
    points(i).period.x=points(i).min.x;
    
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
end

%compute the extra points and features per period
[F,Fdist,points]=compute_features(t,X,points,platThresh,DX);

fcns.compute_features=@compute_features;
fcns.plot_data=@plot_data; %simple plot fcn for one trace of interest
fcns.plot_interactive=@plot_interactive; %adds keypressfcn to switch traces

%plot to show performance
if nargout==0||doPlot==1
    if dokeypress
        plot_interactive(figID,t,X,points)
    else
        if isempty(figID)
            gcf; %new figure if none available, otherwise use current fig
        else
            figure(figID); %bring figID into focus
        end
        plot_data(t,X,points,1)
    end
end

end

function [F,Fdist,points]=compute_features(t,X,points,platThresh,DX)

if ~exist('DX','var')
    DX=slopeY(t,X);
end

nX=size(X,2);

Fdist=repmat( struct('period',0,'apd',0,'spd',0,'pf',0,'amp',[],...
    'baseline',[],'peaks',[],'maxslope',0,'minslope',0,'pthresh',[]) ,1,nX); %'range',[],'thrUp',[],'thrDown',[]

for i=1:nX
    
    nPer=length(points(i).period.t)-1;
    
    if nPer>=1
        
        points(i).period.x=X(points(i).period.ix,i)';
        
        %simple method: baseline=average of successive minima
        baseline=0.5*(points(i).period.x(1:end-1)+points(i).period.x(2:end));
        amp=points(i).max.x-baseline;
        
        pthresh=zeros(1,nPer);
        for j=1:nPer
            ix=points(i).period.ix(j):points(i).period.ix(j+1);
            tt=t(ix);
            xx=X(ix,i);
            dx=DX(ix,i);
            
            thr=baseline(j)+platThresh*amp(j);
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
            
            pthresh(j)=thr;
            
            %TODO: option to recompute max/min (if not same trace as originally used to define periods)
%             [xmax,imax]=max(xx);
%             points(i).max.ix(j)=ix(1)+imax-1;
%             points(i).max.t(j)=tt(imax);
%             points(i).max.x(j)=xmax;
%             
%             [xmin,imin]=min(xx);
%             points(i).min.ix(j)=ix(1)+imin-1;
%             points(i).min.t(j)=tt(imin);
%             points(i).min.x(j)=xmin;
            
            [dxmax,idxmax]=max(dx);
            points(i).dxmax.ix(j)=ix(1)+idxmax-1;
            points(i).dxmax.t(j)=tt(idxmax);
            points(i).dxmax.x(j)=xx(idxmax);
            points(i).dxmax.dx(j)=dxmax;
            
            [dxmin,idxmin]=min(dx);
            points(i).dxmin.ix(j)=ix(1)+idxmin-1;
            points(i).dxmin.t(j)=tt(idxmin);
            points(i).dxmin.x(j)=xx(idxmin);
            points(i).dxmin.dx(j)=dxmin;
        end
        
        
        Fdist(i).period=diff(points(i).period.t); %period defined from minimum to minimum
        Fdist(i).apd=points(i).down.t-points(i).up.t;
        Fdist(i).spd=Fdist(i).period-Fdist(i).apd;
        Fdist(i).pf=Fdist(i).apd./Fdist(i).period;
        Fdist(i).amp=amp;
        Fdist(i).baseline=baseline;
        Fdist(i).peaks=points(i).max.x;
        Fdist(i).maxslope=points(i).dxmax.dx;
        Fdist(i).minslope=points(i).dxmin.dx;
        Fdist(i).pthresh=pthresh;
    else
        %had less than two minima: can't compute features.
        if ~isempty(points(i).min.x)
            Fdist(i).baseline=points(i).min.x;
        end
        if ~isempty(points(i).max.x)
            Fdist(i).peaks=points(i).max.x;
        end
        if length(points(i).max.x)==1 && length(points(i).min.x)==1
            Fdist(i).amp=points(i).max.x-points(i).min.x;
        end
        if ~isempty(points(i).dxmax.dx)
            Fdist(i).maxslope=points(i).dxmax.dx;
        end
        if ~isempty(points(i).dxmin.dx)
            Fdist(i).minslope=points(i).dxmin.dx;
        end
    end
end

%compute summary statistics of feature distributions
fnames=fieldnames(Fdist);

for i=1:length(fnames)
    F.(fnames{i}).mean=mean(Fdist.(fnames{i}));
    F.(fnames{i}).stdev=std(Fdist.(fnames{i}));
end

end


function plot_data(t,X,points,tix)

plot(t,X(:,tix),'k-')
hold on

plot(points(tix).up.t,points(tix).up.x,'bs')
plot(points(tix).down.t,points(tix).down.x,'bo')
plot(points(tix).max.t,points(tix).max.x,'rv')
plot(points(tix).min.t,points(tix).min.x,'r^')
plot(points(tix).dxmax.t,points(tix).dxmax.x,'g>')
plot(points(tix).dxmin.t,points(tix).dxmin.x,'g<')
% plot(xlim(),thrUp(tix)*[1,1],'r--')
% if thrUp~=thrDown
%     plot(xlim(),thrDown(tix)*[1,1],'b--')
% end
plot(points(tix).period.t,points(tix).period.x,'bs')
hold off

xlabel('t')
ylabel('x')
axis tight
YLIM=ylim();
ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
end

function plot_interactive(figID,t,X,points)
tix=1;
if isempty(figID)
    figID=gcf; %new figure if none available, otherwise current fig
else
    figID=figure(figID);
end
figID.KeyPressFcn=@keypressFcn;
figID.UserData.tix=tix; %store tix (trace to plot) in userdata
figID.UserData.t=t;
figID.UserData.X=X;
figID.UserData.points=points;

plot_data(t,X,points,tix)
end

function keypressFcn(src,event)
tix=src.UserData.tix;
t=src.UserData.t;
X=src.UserData.X; nX=size(X,2);
points=src.UserData.points;
switch(event.Key)
    case {'leftarrow'}
        if tix>1
            tix=tix-1;
        end
    case {'rightarrow'}
        if tix<nX
            tix=tix+1;
        end
end
plot_data(t,X,points,tix)
src.UserData.tix=tix;
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