classdef Experiment < handle
    %class to represent a single experiment containing >=1 timeseries that share the same time sample points. For
    %example, a fluorescence videomicroscopy experiment with multiple ROIs in a field of view. 

    %TODO: plotting efficiently - plot all timeseries onces, and use an
    %update function to set the visibility on/off (also would support
    %option to plot out-of-focus traces as light gray)
    
    properties
        name=''
        date=''
        sex=''
        condition='' %eg. WT vs KO
        filename=''
        notes={}
        
        t
        dt
        
        %segment is a subinterval of the trace, the basic unit within which
        %periodicity will be detected and features will be computed
        segment=struct('name','','endpoints',[],'ix',[],...
            'points',struct(),'features_periods',struct(),'features_trace',struct())
        nS
        
        %TODO: support 3D array - 3rd dim is observable id (one [nT x nX] page per observable)
        % this requires supporting different pipeline options for each
        % observable and allows computation of new higher order features,
        % such as phase relationship
        % --> actually, will likely determine period markers using only one
        % observable, then compute features relative to those markers in
        % other traces?
        
        %TODO: Customizeable pipeline? eg. array processingStep objects to define
        %steps (name+operation+method+params), and pointer to which one
        %gets periods measured and which one features get measured. pages
        %of X correspond to each step: eg. X(:,:,1)=raw, X(:,:,2)=norm, ..., X(:,:,end)=per
        X
        Xnorm
        Xtrend
        Xdetrend
        Xfilt
        
        include %set element to zero to exclude a trace; X(:,include)
        nX
        
        featureFcn
        fnames_periods
        fnames_trace
        
        fs %1/dt
        f
        psd
        
        group
        nG
        Xg %average traces with same group (always use Xdetrend? probably)
        Xgfilt %filter using filterMethod, filterParam
    end
    
    %store parameters used
    properties
        
%         tOmit
%         nOmit

%support variable length parameter lists, eg. for name-value pairs
        normMethod='none'
        normParam=[]
        trendMethod='none'
        trendParam=[]
        filterMethod='none'
        filterParam=[]
        
        averageMethod='arithmetic'
        averageParam=[]
        
        featureMethod='threshold'
        featureParam={}
        
        interpMethod='linear'
        
%         thrFrac=[0.55,0.45]
%         delta=0.25
%         Tbig=15
%         Tsmall=4
        
        %keep track of which trace is in focus for plotting - no, this should be a property of the figure
%         tix=1;
    end
    
    methods
        %constructor
        % Experiment(): requests file
        % Experiment(filename): looks for given file
        % Experiment(data): uses numeric matrix data=[t,X]
        % TODO: Experiment(___,Name,Value) : name-value pairs to set parameters, metadata, etc.
        function expt=Experiment(varargin)
            
            %parse inputs
            loadFile=true;
            fullfile=[];
            if ~isempty(varargin) && ( ischar(varargin{1}) || isstring(varargin{1}) )
                [path,filename,ext]=fileparts(varargin{1});
                fullfile=[path,filesep,filename,ext];
            elseif ~isempty(varargin) && isnumeric(varargin{1})
                data=varargin{1};
                loadFile=false;
            end
            
            %load a new file
            if loadFile
                if isempty(fullfile) || ~exist(fullfile,'file')
                    [filename,path]=uigetfile({'*.xls*'});
                    fullfile=[path,filename];
                end
                
                [num,txt,raw]=xlsread(fullfile);
                %extract data
                ixs=find(cellfun(@(x)(isnumeric(x)&&~isnan(x)),raw(:,1)),1,'first');
                data=cell2mat(raw(ixs:end,1:size(num,2)));
                
                %look for special metadata entries
                headerNames=txt(:,1);
                keyword={'name','date','sex','condition'};
                for i=1:length(keyword)
                    r=find(strcmpi(headerNames,keyword{i}));
                    if ~isempty(r)
                        expt.(keyword{i})=txt{r,2};
                    end
                end
                
                r=find(strcmpi(headerNames,'group'));
                if ~isempty(r)
                    expt.setGroup(cell2mat(raw(r,2:size(num,2))))
                end
                
                expt.filename=fullfile;
            end
            
            time=data(:,1);
            time=time-time(1); %shift time to start at zero
            Xraw=data(:,2:end);
            
            DT=mode(diff(time)); 

            %interpolate any missing data using neighboring time points
            %TODO: make this optional?
            
            timeIsNan=isnan(time); %this would be a problem!
            if any(timeIsNan)
                expt.notes=[expt.notes,{'missing value in time column'}];
            end
            
            skip=any(abs(diff(time)-DT)>2*DT);
            if any(skip)
                expt.notes=[expt.notes,{'large timestep detected'}];
            end
            
            rowIsNan=any(isnan(Xraw),2);
            if any(rowIsNan)
                expt.notes=[expt.notes,{'some rows have NaN'}];
                r=find(rowIsNan);

                if r(1)==1
                    time=time(2:end,1);
                    Xraw=Xraw(2:end,:);
                    r=r(2:end)-1;
                end
                for ii=1:length(r)
                    for jj=1:size(Xraw,2)
                        iix=[r(ii)-1,r(ii)+1];
                        x=time(iix,1);
                        v=Xraw(iix,jj);
                        xq=time(r(ii),1);
                        Xraw(r(ii),jj)=interp1(x,v,xq,expt.interpMethod);
                    end
                end
            end
%             expt.t=time;
%             expt.X=Xraw;
            
            %interpolate on a fixed grid with given DT
            %TODO: DT as parameter (for same DT across expts)
            % - alt, make this a method with DT parameter
            XI=[];
            ti=0:DT:time(end); ti=ti';
            for j=1:size(Xraw,2)
                XI(:,j)=interp1(time,Xraw(:,j),ti,'pchip');
            end
            expt.t=ti;
            expt.X=XI;
            
            expt.nX=size(expt.X,2);
            
            expt.dt=DT;
            expt.fs=1/expt.dt;
            
            %default segment = whole trace
            expt.defineSegments({'1'},[ti(1),ti(end)]);
%             expt.setGroup(ones(1,expt.nX));
            expt.include=true(1,expt.nX);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup helper functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function defineSegments(expt,names,startTimes)
            %expects names as cell array, times as vector
            expt.nS=length(names);
%             expt.segment(expt.nS).name=''; %expand the array of struct
            startTimes(end+1)=expt.t(end); %append final time
            for i=1:expt.nS
                expt.segment(i).name=names{i};
                expt.segment(i).endpoints=[startTimes(i),startTimes(i+1)];
                expt.segment(i).ix=expt.t>=startTimes(i) & expt.t<=startTimes(i+1);
            end
        end
        
        function setGroup(expt,groupvar)
            expt.group=groupvar;
            expt.nG=unique(groupvar);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preprocessing functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TODO: forced order, but flags to choose which steps (if any) to apply?
        %     - alt: arbitrary list of operations, with fully qualified options.
        %TODO: per-interval, or whole trace options (eg. linear detrend per interval)
        %TODO: single preprocessing convenience function. how to handle all the options+params? 
        
        function normalize(expt,method,methodPar,doPlot)
            expt.Xnorm=normalizeTraces(expt.t,expt.X,method,methodPar,doPlot);
            expt.normMethod=method;
            expt.normParam=methodPar;
        end
        
        %detrend supports per-interval: for piecewise linear detrend.
        function detrend(expt,method,methodPar,perSegment,doPlot)
            if isempty(expt.Xnorm)
                expt.Xnorm=expt.X;
                expt.normMethod='none';
                expt.normParam=[];
            end
            if perSegment
                for i=1:expt.nS
                    thisIx=expt.segment(i).ix;
                    %TODO: this plots a figure for each segment. could suppress & make combined plot?
                    [Xdt,Xt]=detrendTraces(expt.t(thisIx),expt.Xnorm(thisIx,:),method,methodPar,doPlot);
                    %TODO: make segments join
                    expt.Xtrend(thisIx,:)=Xt;
                    expt.Xdetrend(thisIx,:)=Xdt;
                end
            else
                [expt.Xdetrend,expt.Xtrend]=detrendTraces(expt.t,expt.Xnorm,method,methodPar,doPlot);
            end
            expt.trendMethod=method;
            expt.trendParam=methodPar;
        end
        
        function averageTraces(expt)
            XDTg=zeros(length(expt.t),expt.nG);
            for j=1:nI  
                XDTg(:,j)=mean(XDTg(:,I==uI(j)),2);  
            end
            expt.Xg=XDTg;
        end
        
        
        function filter(expt,method,methodPar,doPlot)
            if isempty(expt.Xnorm)
                expt.Xnorm=expt.X;
                expt.normMethod='none';
                expt.normParam=[];
            end
            if isempty(expt.Xdetrend)
                expt.Xdetrend=expt.X;
                expt.trendMethod='none';
                expt.trendParam=[];
            end
            expt.Xfilt=filterTraces(expt.t,expt.Xdetrend,method,methodPar,doPlot);
            expt.filterMethod=method;
            expt.filterParam=methodPar;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TODO expand this - detect_periods return only points/fcns, then compute_features should apply fcn to relevant
        %traces.
        function compute_features(expt,method,varargin)
            expt.detect_periods(method,varargin{:});
            expt.periodogram();
            expt.fnames_periods=fieldnames(expt.segment(1).features_periods)';
            expt.fnames_trace=fieldnames(expt.segment(1).features_trace)';
            expt.featureMethod=method;
            expt.featureParam=varargin; %methodpar should be varargin?
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plotTrace(expt,whichPlot,tix,showPts)
            %dispatch function for plotting expt data
                        
            if ~exist('showPts','var')||isempty(showPts)
                showPts=true;
            end
            
            %no tix => interactive figure
            if ~exist('tix','var')||isempty(tix)
                tix=1;
                figID=gcf; %newfig?
                figID.KeyPressFcn={@expt.traceKeypress,whichPlot};
                figID.UserData=tix; %store the current trace ID with figure
            end
            
            expt.plot_t(whichPlot,tix,showPts)
            
        end
        
        function plotPeriodogram(expt,tix)
            %overlay psd for each segment
            if ~isempty(expt.psd)
                       
            if ~exist('tix','var')||isempty(tix)
                tix=1;
                figID=gcf; %newfig?
                figID.KeyPressFcn=@expt.psdKeypress;
                figID.UserData=tix; %store the current trace ID with figure
            end
            
            expt.plot_psd(tix);
            
            end
        end
        
        function hl=plotFeatures(expt,featureType,xname,yname,tix)
            %TODO: how to manage plotting options?
                
%             if ~isempty(expt.segment(1).features)
                            
            if isempty(xname)
                xname='segment';
            end
            
            if ~exist('tix','var')||isempty(tix)
                tix=1;
                figID=gcf; %newfig?
                figID.KeyPressFcn={@expt.featureKeypress,xname,yname};
                figID.UserData=tix; %store the current trace ID with figure
            end
            
            switch featureType
                case 'periods'
                    xname=validatestring(xname,['t','segment',expt.fnames_periods]);
                    yname=validatestring(yname,expt.fnames_periods);
                    hl=expt.plot_features_periods(xname,yname,tix);
                case 'trace'
                    xname=validatestring(xname,['t','segment',expt.fnames_trace]);
                    yname=validatestring(yname,expt.fnames_trace);
                    hl=expt.plot_features_trace(xname,yname,tix);
            end
            
            axis tight
            YLIM=ylim();
            ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
            
            switch xname
                    case {'t'}
                        xlabel('t')
                        for i=1:expt.nS
                            if strcmpi(class(hl(i)),'matlab.graphics.chart.primitive.Line')
                            hl(i).LineStyle='-';
                            end
                        end
                        hold on
                        for i=2:expt.nS
                            plot(expt.segment(i).endpoints(1)*[1,1],ylim(),'g')
                        end
                        hold off
                        xlim([expt.t(1),expt.t(end)])
                        
                    case {'segment'} %TODO: do jitter here?
                        xticks(1:expt.nS)
                        xticklabels({expt.segment.name})
                        xlim([0.5,expt.nS+0.5])
                        xlabel('segment')
                        
                otherwise
                        xlabel(xname)
            end
            
            ylabel(yname)
            
%             end
        end
        
        %save results. If required metadata is not set, force user to enter it now
        function save_results(expt)
        end
        
    end
    
    methods (Access=private)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function detect_periods(expt,method,varargin)
            
            extras={};
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                tt=expt.t(thisIx);
                xx=expt.Xfilt(thisIx,:);
                switch method
                    case {'peaks'}
                        delta=varargin{1};
                        if length(varargin)>1, extras=varargin(2:end); end
                        [F, Fdist, points, fcns]=peak_detector(tt,xx,delta,extras{:});

                    case {'threshold'}
                        frac=varargin{1};
                        if length(varargin)>1, extras=varargin(2:end); end
                        [F, Fdist, points, fcns]=plateau_detector(tt,xx,frac,extras{:});
                end
                
                expt.segment(i).points=points;
                expt.segment(i).features_periods=Fdist;
                expt.segment(i).features_trace=F;
            end
            expt.featureFcn=fcns.compute_features;
        end
        
        function periodogram(expt)
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                [PSD,F,Pmax,fmax]=powerSpectrum(expt.Xdetrend(thisIx,:),expt.fs);
%                 [PSD,F,Pmax,fmax]=powerSpectrum(expt.Xfilt(thisIx,:),expt.fs);
                Tpsd=1./fmax;
                expt.psd(:,:,i)=PSD;
                
                Tpsd=num2cell(Tpsd);
                fmax=num2cell(fmax);
                Pmax=num2cell(Pmax);
                [expt.segment(i).features_trace.Tpsd]=Tpsd{:};
                [expt.segment(i).features_trace.fmax]=fmax{:};
                [expt.segment(i).features_trace.Pmax]=Pmax{:};
            end
            expt.f=F; %frquency vector
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot_t(expt,whichPlot,tix,showPts)
            
            %BUG: matlab errors when mousing over data. datatips?? super
            %annoying
            
            %TODO: support a line handle passed in instead of ax? i.e. just
            %change x/y data of the line. still need to remove old pts and
            %plot new ones if desired.
            
%             axes(ax)
            cla
            hold on
            
            points_option='period';
            x2=[];
            switch whichPlot
                case {'raw'}
                    x=expt.X(:,tix);
                    
                case {'norm'}
                    x=expt.Xnorm(:,tix);
                    x2=expt.Xtrend(:,tix);
                    
                case {'detrend'}
                    x=expt.Xdetrend(:,tix);
                    x2=expt.Xfilt(:,tix);
                    
                case {'filt'}
                    x=expt.Xfilt(:,tix);
                    points_option='all';
                    
                otherwise
                    error([whichPlot, ' is not a supported trace to plot']);
            end
            
            plot(expt.t,x,'k')
            if ~isempty(x2)
                plot(expt.t,x2)
            end
            
            
            %slower:
%             if ~showPts
%                 points_option='none';
%             end
%             plotTrace=false;
%             for i=1:expt.nS
%                 tt=expt.t(expt.segment(i).ix);
%                 xx=x(expt.segment(i).ix);
%                 pts=expt.segment(i).points(tix);
%                 pts.period.x=interp1(tt,xx,pts.period.t); %interp x value of period markers
%                 
%                 expt.segment(i).plot(tt,xx,pts,1,plotTrace,points_option);
%             end

            if showPts
                for i=1:expt.nS
                    tt=expt.t(expt.segment(i).ix);
                    xx=x(expt.segment(i).ix);
                    pts=expt.segment(i).points(tix);
                    
                    if length(pts.period.t)>1
                        
                    if points_option=="all"
                        plot(pts.min.t,pts.min.x,'r^')
                        plot(pts.max.t,pts.max.x,'rv')
                        plot(pts.dxmin.t,pts.dxmin.x,'g<')
                        plot(pts.dxmax.t,pts.dxmax.x,'g>')
                        plot(pts.up.t,pts.up.x,'bd')
                        plot(pts.down.t,pts.down.x,'bo')
%                         line(pts.min.t,pts.min.x,'color','r','marker','^','linestyle','none')
%                         line(pts.max.t,pts.max.x,'color','r','marker','v','linestyle','none')
%                         line(pts.dxmin.t,pts.dxmin.x,'color','g','marker','<','linestyle','none')
%                         line(pts.dxmax.t,pts.dxmax.x,'color','g','marker','>','linestyle','none')
%                         line(pts.up.t,pts.up.x,'color','b','marker','d','linestyle','none')
%                         line(pts.down.t,pts.down.x,'color','b','marker','o','linestyle','none')
                        
%                         tupdwn=[pts.up.t(1:length(pts.down.t)); pts.down.t];
%                         ythresh=[1;1]*expt.segment(i).features(tix).pthresh;
%                         plot(ax,tupdwn,ythresh,'b-')
                    else
                    
                        xPer=interp1(tt,xx,pts.period.t); %this is necessary except for filt
                        plot(pts.period.t,xPer,'bs')
%                         line(pts.period.t,xPer,'color','b','marker','s','linestyle','none')

                    end
                    
                    end
                end
            end
            
            axis tight
%             xlabel('t')
%             ylabel(['X',whichPlot])
            YLIM=ylim();
            ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
            
            for i=2:expt.nS
                plot(expt.segment(i).endpoints(1)*[1,1],ylim(),'g')
            end
            
            hold off
        end
        
        function plot_psd(expt,tix)
            
%             axes(ax)
            cla
            hold on;
            
            for i=1:expt.nS
                plot(expt.f,expt.psd(:,tix,i)); 
%                 plot(expt.f,pow2db(expt.psd(:,tix,i)));
            end
            % set(gca,'yscale','log')
%             fhi=find(
%             title('One Sided Power Spectral Density');       
            xlabel('frequency')         
            ylabel('power');
            axis tight
%             xlim(ax,[0,0.5])
            if expt.nS>1
                legend({expt.segment.name})
            end
            hold off
        end
        
        function hl=plot_features_periods(expt,xname,yname,tix)
            
            cla
            hold on;
            
            xJitter=0.05; %make this a property of the class?
           
            XX=cell(1,expt.nS);
            YY=cell(1,expt.nS);
            colorBySegment=false;  %TODO: always color by segment?
            for i=1:expt.nS
                        
                y=expt.segment(i).features_periods(tix).(yname);
                if isempty(y)
                    return
                end
                
                switch xname
                    case {'t'}
                        %TODO: option for marker location?
                         x=(expt.segment(i).points(tix).period.t(1:end-1)+expt.segment(i).points(tix).period.t(2:end))/2; %midpoint of period
%                         x=expt.segment(i).points(tix).max.t;
%                         x=expt.segment(i).points(tix).min.t(1:end-1); %first min
%                         x=expt.segment(i).points(tix).min.t(2:end); %second min
%                         x=expt.segment(i).points(tix).down.t;
%                         x=expt.segment(i).points(tix).up.t;
%                         x=(expt.segment(i).points(tix).down.t+expt.segment(i).points(tix).period.t(2:end))/2; %midpoint of silent post active
%                         x=(expt.segment(i).points(tix).up.t+expt.segment(i).points(tix).period.t(1:end-1))/2; %midpoint of silent pre active   
                        
                    case {'segment'}
                        %TODO: option for violin plot, boxplot?
                        x=i+xJitter*randn(1,length(y));

                    case expt.fnames
                        x=expt.segment(i).features(tix).(xname);
                        colorBySegment=true;

                    otherwise
                        error(['Invalid name for x-axis: ' xname])
                end
                
                if isempty(x)
                    y=[];
                end
                
                XX{i}=x;
                YY{i}=y;

            end
                          
            
            ax=gca;
            ax.ColorOrderIndex=1;  
            hl=matlab.graphics.chart.primitive.Line.empty(expt.nS,0);
            for i=1:expt.nS
                if ~isempty(XX{i})
                    hl(i,1)=plot(XX{i},YY{i},'o');
                end
            end
            if colorBySegment
                legend(hl,{expt.segment.name}) %TODO: option to suppress legend?
            else
                for i=1:expt.nS
                    if ~isempty(XX{i})
                        hl(i).Color='k';
                    end
                end
            end
            
            hold off
            
        end
        
        function hl=plot_features_trace(expt,xname,yname,tix)
            
%             axes(ax)
            cla
            hold on;
            
            xJitter=0.05;
           
            HX=[];
            HY=[];
            XX=[];
            YY=[];
            
            colorBySegment=false;
            for i=1:expt.nS
                
                switch xname
                        
                    case {'segment'}
                        %plot all trace data as line seg, highlight point tix
%                         xx=i+xJitter*randn(1,expt.nX);
                        xx=i*ones(1,expt.nX);
                        yy=[expt.segment(i).features_trace.(yname)];
                        HX(end+1,1)=i;
                        HY(end+1,1)=yy(tix);
                        

                    case expt.fnames_trace
                        xx=[expt.segment(i).features_trace.(xname)];
                        yy=[expt.segment(i).features_trace.(yname)];
                        colorBySegment=true;
                        HX(end+1,1)=xx(tix);
                        HY(end+1,1)=yy(tix);

                    otherwise
                        error(['Invalid name for x-axis: ' xname])
                end
                
                XX(end+1,:)=xx;
                YY(end+1,:)=yy;

            end
            
            %line segment below markers:
            hlAll=[];
            hlTix=[];
            if expt.nS>1
                hlAll=plot(XX,YY,'-','color',[0.5,0.5,0.5]);
                hlTix=plot(HX,HY,'k-','linewidth',1.5);
            end
            
            %markers
            
            ax=gca;
            ax.ColorOrderIndex=1;
            hmAll=plot(XX',YY','o');
            for i=1:expt.nS
                hmTix(i)=plot(HX(i),HY(i),'o');
                hmTix(i).Color=hmAll(i).Color;
                hmTix(i).MarkerFaceColor=hmAll(i).Color;
            end
            if expt.nS>1
                legend(hmAll,{expt.segment.name})
            end
            
            hold off
            
            hl.hlAll=hlAll;
            hl.hlTix=hlTix;
            hl.hmAll=hmAll;
            hl.hmTix=hmTix;
        end
        
        function traceKeypress(expt,src,event,whichPlot,showPts)
            tix=src.UserData;
            switch(event.Key)
                case {'leftarrow'}
                    if tix>1
                        tix=tix-1;
                    end
                case {'rightarrow'}
                    if tix<expt.nX
                        tix=tix+1;
                    end
            end
            expt.plot_t(whichPlot,tix,showPts)
            src.UserData=tix;
        end
        
        function psdKeypress(expt,src,event)
            tix=src.UserData;
            switch(event.Key)
                case {'leftarrow'}
                    if tix>1
                        tix=tix-1;
                    end
                case {'rightarrow'}
                    if tix<expt.nX
                        tix=tix+1;
                    end
            end
            expt.plot_psd(tix)
            src.UserData=tix;
        end
        
        function featureKeypress(expt,src,event,xname,yname)
            tix=src.UserData;
            switch(event.Key)
                case {'leftarrow'}
                    if tix>1
                        tix=tix-1;
                    end
                case {'rightarrow'}
                    if tix<expt.nX
                        tix=tix+1;
                    end
            end
            expt.plot_features(xname,yname,tix);
            src.UserData=tix;
        end
    end

end