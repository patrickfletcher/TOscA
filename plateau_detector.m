function [points,features]=plateau_detector(t, X, varargin)
%PLATEAU_DETECTOR A plateau detector for oscillating timeseries. Uses Shmitt
% trigger concept - upward transition threshold >= downward transition
% threshold. The timeseries is divided into periods measured from upward
% transition to upward transition. 
%
% If results are plotted, use left/right arrows to switch between columns
% of X.
%
% Required Inputs:
%  t - vector of timepoints
%  X - vector matching t, or matrix with columns matching t
%
% Optional input:
%  f - fraction of amplitude for threshold: 
%     scalar {default, f=0.5}, upward=downward threshold,
%     vector, f=[upFraction,downFraction] with upFraction>=downFraction.
%
% optional name-value pair arguments:
%  'Interpolate', {true}/false
%  'ThresholdPercentiles', default=[0,100] (threshold is fraction of this range of amplitudes)
%  'Plot',true/{false} (use to control plotting with output arguments)
%
% Outputs:
%  features - struct containg feature distributions for each column of X
%  points - struct containing special points, for plotting

%check inputs and parse optional parameters

periodMarker='min';

if isvector(X)
    X=X(:);
end
    
[thrPtiles,fracUp,fracDown,minAmp,threshMethod,doInterp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

nX=size(X,2); %number of traces

globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

sufficientGlobalAmp=globalXamp>minAmp;

%independent thresholds for each trace
switch lower(threshMethod)
    case 'fracamp'
        %fraction of global amplitude method
        thrUp=globalXmin + fracUp*globalXamp;
        if fracDown>0
            thrDown=globalXmin + fracDown*globalXamp;
        else
            thrDown=thrUp;
        end
        
    case 'rawvals'
        %specific raw values as thresholds
        thrUp=fracUp;
        thrDown=fracDown;
        
    case 'medfrac'
        %median +/- fraction of amp
        medX=median(X,1);
        thrUp=medX + fracUp*globalXamp;
        thrDown=medX - fracDown*globalXamp;
        
    case 'medraw'
        %median +/- raw values
        medX=median(X,1);
        thrUp=medX + fracUp;
        thrDown=medX - fracDown;
end


%TODO: other possibilities
% central value +/- fraction of global amp (eg. median +/- 0.05*amp)

%option for dimensional values as thresholds for above methods


%TODO: make sure initialization of upstate toggle is robust
DX=slopeY(t,X); %todo: noise-robust method; measure between tmin(i) and tmin(i+1)?
% isUp=ones(1,nX); %forces first point to be downward transition
% isUp=zeros(1,nX); %forces first point to be upward transition
isUp=X(1,:)>thrUp;
% isUp=X(1,:)>=thrUp | (X(1,:)<=thrUp&X(1,:)>=thrDown&DX(1,:)<0);

pt=struct('ix',[],'t',[],'x',[],'dx',[]);

points=repmat(struct('period',pt,'up',pt,'down',pt,...
    'max',pt,'min',pt,'dxmax',pt,'dxmin',pt),1,nX);


for i=2:length(t)
    
    %if strated up, then went below thrUp but not below thrDown, set isUp=0
%     firstDown=isUp & X(i,:)<thrUp & arrayfun(@(x)isempty(x.up.ix),points);
%     isUp(firstDown)=0;
    
    %check for downward-transition and upward-transition
    isDownTransition=isUp&X(i,:)<=thrDown & sufficientGlobalAmp;
    isUpTransition=~isUp&X(i,:)>=thrUp & sufficientGlobalAmp;
    
    if any(isDownTransition)
        idx=find(isDownTransition);
        for j=idx
            %index
            points(j).down.ix(end+1)=i;
            
            if doInterp
                %interpolate t
                thisX=thrDown(j);
                thisT=(thrDown(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            points(j).down.t(end+1)=thisT;
            points(j).down.x(end+1)=thisX;
        end
        isUp(isDownTransition)=0;
    end
    
    if any(isUpTransition)
        idx=find(isUpTransition);
        for j=idx
            points(j).up.ix(end+1)=i;
            
            if doInterp
                %interpolate t
                thisX=thrUp(j);
                thisT=(thrUp(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            points(j).up.t(end+1)=thisT;
            points(j).up.x(end+1)=thisX;
        end
        isUp(isUpTransition)=1;
    end
    
end

numUp=arrayfun(@(x) length(x.up.t),points);

features=struct('period',[],'APD',[],'PF',[],'amp',[],'baseline',[],'peaks',[],'pthresh',[],'maxslope',[],'minslope',[]);
features=struct();
for i=1:nX
        
    if numel(points(i).up.t)>1
    
%     trim
    if points(i).down.t(1)<=points(i).up.t(1)
        points(i).down.ix=points(i).down.ix(2:end);
        points(i).down.t=points(i).down.t(2:end);
        points(i).down.x=points(i).down.x(2:end);
    end
    if points(i).down.t(end)>=points(i).up.t(end)
        points(i).down.ix=points(i).down.ix(1:end-1);
        points(i).down.t=points(i).down.t(1:end-1);
        points(i).down.x=points(i).down.x(1:end-1);
    end
    
    points(i).period.ix=points(i).up.ix;
    points(i).period.t=points(i).up.t;
    points(i).period.x=points(i).up.x;
    
    nT=length(points(i).up.t)-1;
    
    features(i).period=diff(points(i).up.t);
    features(i).APD=points(i).down.t-points(i).up.t(1:end-1); %active phase duration
    features(i).PF=features(i).APD./features(i).period;
    
    for j=1:nT
        tt=points(i).up.ix(j):points(i).up.ix(j+1)-1;
        [xmax,imax]=max(X(tt,i));
        points(i).max.ix(j)=tt(1)+imax-1;
        points(i).max.t(j)=t(tt(imax));
        points(i).max.x(j)=xmax;
        
        [xmin,imin]=min(X(tt,i));
        points(i).min.ix(j)=tt(1)+imin-1;
        points(i).min.t(j)=t(tt(imin));
        points(i).min.x(j)=xmin;
        
        %compute min slope anywhere in the period
        [dxmin,idxmin]=min(DX(tt,i));
        points(i).dxmin.ix(j)=tt(1)+idxmin-1;
        points(i).dxmin.t(j)=t(tt(idxmin));
        points(i).dxmin.x(j)=X(tt(idxmin),i);
        points(i).dxmin.dx(j)=dxmin;
        
        %compute maxslope only in the active phase? 
        tt=points(i).up.ix(j):points(i).down.ix(j);
        [dxmax,idxmax]=max(DX(tt,i));
        points(i).dxmax.ix(j)=tt(1)+idxmax-1;
        points(i).dxmax.t(j)=t(tt(idxmax));
        points(i).dxmax.x(j)=X(tt(idxmax),i);
        points(i).dxmax.dx(j)=dxmax;
    end
    
    features(i).baseline=points(i).min.x;
    features(i).peaks=points(i).max.x;
    features(i).amp=points(i).max.x-points(i).min.x;
    features(i).range=globalXamp(i);
    features(i).thrUp=thrUp(i);
    features(i).thrDown=thrDown(i);
    features(i).pthresh=(thrUp(i)+thrDown(i))/2;
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
    if isempty(figID)
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

        plot(t,X(:,tix),'k-','Tag','plateau_detector')
        axis tight
        hold on
        plot(points(tix).up.t,points(tix).up.x,'bs','Tag','plateau_detector')
        plot(points(tix).down.t,points(tix).down.x,'bo','Tag','plateau_detector')
        plot(points(tix).max.t,points(tix).max.x,'rv','Tag','plateau_detector')
        plot(points(tix).min.t,points(tix).min.x,'r^','Tag','plateau_detector')
        plot(points(tix).dxmax.t,points(tix).dxmax.x,'g>','Tag','plateau_detector')
        plot(points(tix).dxmin.t,points(tix).dxmin.x,'g<','Tag','plateau_detector')
        plot(xlim(),thrUp(tix)*[1,1],'r--','Tag','plateau_detector')
        if thrUp~=thrDown
            plot(xlim(),thrDown(tix)*[1,1],'b--','Tag','plateau_detector')
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

function [thrPtiles,fracUp,fracDown,minAmp,threshMethod,doInterp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin)

    %default parameters
    defaultF=[0.5,0.4]; %near halfmax
    defaultthrPtiles=[0,100];
    defaultminAmp=0;
    defaultThresholdMethod='medfrac';
    doInterp=true;
    doPlot=false;
    doKeypress=true;
    figID=[];
    
    p=inputParser;
    validX=@(x) isreal(x) && size(x,1)==length(t); %traces are columns of X
    validSwitch=@(x) isscalar(x) && (isnumeric(x)||islogical(x));
    addRequired(p,'t',@(x) isreal(x));
    addRequired(p,'X',validX);
    addOptional(p,'f',defaultF);  %somehow adding the validation function here messes things up
    addParameter(p,'ThresholdPercentiles',defaultthrPtiles,@(x) isreal(x) && numel(x)==2);
    addParameter(p,'MinimumAmplitude',defaultminAmp,@isreal);
    addParameter(p,'ThresholdMethod',defaultThresholdMethod);
    addParameter(p,'Interpolate',doInterp,validSwitch);
    addParameter(p,'Plot',doPlot,validSwitch);
    addParameter(p,'Keypress',doKeypress,validSwitch);
    addParameter(p,'FigureID',figID);
    
    parse(p,t,X,varargin{:});
    
    f=p.Results.f;
    thrPtiles=p.Results.ThresholdPercentiles;
    minAmp=p.Results.MinimumAmplitude;
    threshMethod=p.Results.ThresholdMethod;
    doInterp=p.Results.Interpolate;
    doPlot=p.Results.Plot;
    figID=p.Results.FigureID;
    dokeypress=p.Results.Keypress;
    
    if isscalar(f)
        fracUp=f;
        fracDown=f;
    elseif length(f)==2
        fracUp=f(1);
        fracDown=f(2);
    else
        error('Invalid value for fractional threshold')
    end
    
%     %snap to fracUp?
%     if fracDown>fracUp
%         fracDown=fracUp;
%     end
    
end