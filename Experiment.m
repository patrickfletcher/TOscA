classdef Experiment < handle
    %class to represent a single experiment containing >=1 timeseries that share the same time sample points. For
    %example, a fluorescence videomicroscopy experiment with multiple ROIs in a field of view. 

    properties
        name=''
        date=''
        sex=''
        animalCondition='' %eg. WT vs KO
        filename=''
        notes={}
        
        t
        dt
        
        segment=struct('name','','endpoints',[],'ix',[])
        nS
        
        %TODO: support 3D array - 3rd dim is observable id (one [nT x nX] page per observable)
        X
        Xnorm
        Xtrend
        Xdetrend
        Xfilt
        
        include %set element to zero to exclude a trace; X(:,include)
        nX
        
        %don't want Expt to have to know
        %if saving undeclared struct doesn't work try cell array for each segment
        %maybe the detector functions could export struct definition when called with no inputs
        points={}
        features={}
        fnames
        
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

        normMethod='none'
        normParam=[]
        trendMethod='none'
        trendParam=[]
        filterMethod='none'
        filterParam=[]
        
        averageMethod='arithmetic'
        averageParam=[]
        
        interpMethod='linear'
        
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
            header={};
            if loadFile
                if isempty(fullfile) || ~exist(fullfile,'file')
                    [filename,path]=uigetfile({'*.xls*'});
                    fullfile=[path,filename];
                end
                
                [data,header]=xlsread(fullfile);
                expt.filename=fullfile;
                expt.name=filename;
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
            expt.defineSegments({''},[ti(1),ti(end)]);

            expt.include=true(1,expt.nX);
        end
        
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
        
        %save results. If required metadata is not set, force user to enter it now
        function save()
        end
        
        %Preprocessing functions
        %TODO: forced order?
        %TODO: per-interval, or whole trace option (eg. linear detrend per interval)
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
        
        
        %Analysis functions - per interval
        function periodogram(expt)
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                [PSD,F,pmax,fmax]=powerSpectrum(expt.Xfilt(thisIx,:),expt.fs);
                Tpsd=1./fmax;
                expt.psd(:,:,i)=PSD;
                
                for j=1:expt.nX
                    expt.features{i}(j).Tpsd=Tpsd(j);
                    expt.features{i}(j).fmax=fmax(j);
                    expt.features{i}(j).Pmax=pmax(j);
                end
            end
            expt.f=F;
        end
        
        function detect_periods(expt,method,methodPar)
%             expt.points={};
%             expt.features={};
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                tt=expt.t(thisIx);
                xx=expt.Xfilt(thisIx,:);
                switch method
                    case {'peaks'}
                        delta=methodPar;
                        [pts,feats]=peak_detector(tt,xx,delta);

                    case {'threshold'}
                        frac=methodPar;
                        [pts,feats]=plateau_detector(tt,xx,frac);
                end
                expt.points{i}=pts;
                expt.features{i}=feats;
                
                %also get xvalues for all X types
                %independently max/min per period for non-filt traces (detrend at least)
                
                %get average features per segment per trace
                
            end
            expt.fnames=fieldnames(feats)';
        end
        
        
        function plotTrace(expt,ax,whichPlot,showPts,tix)
            %dispatch function for plotting expt data
                        
            if isempty(ax)
                ax=gca;
            end
            
            if ~exist('showPts','var')||isempty(showPts)
                if ~isempty(expt.points)
                    showPts=true;
                else
                    showPts=false;
                end
            end
            
            %no tix => interactive figure
            if ~exist('tix','var')||isempty(tix)
                tix=1;
                figID=gcf; %newfig?
                figID.KeyPressFcn={@expt.traceKeypress,whichPlot,showPts};
                figID.UserData=tix; %store the current trace ID with figure
                ax=gca;
            end
            
            expt.plot_t(ax,whichPlot,tix,showPts)
            
        end
        
        function plotPeriodogram(expt,ax,tix)
            %overlay psd for each segment
            if isempty(expt.psd)
                expt.periodogram();
            end
                       
            if isempty(ax)
                ax=gca;
            end
            
            if ~exist('tix','var')||isempty(tix)
                tix=1;
                figID=gcf; %newfig?
                figID.KeyPressFcn=@expt.psdKeypress;
                figID.UserData=tix; %store the current trace ID with figure
                ax=gca;
            end
            
            expt.plot_psd(ax,tix);
            
        end
        
        function plotFeature(expt,ax,xname,yname,tix)
            %xname={t,segment,fnames}
            %yname={fnames}
            
        end
        
    end
    
    methods (Access=private)
        
        function plot_t(expt,ax,whichPlot,tix,showPts)
            
            doAllPts=false;
            switch whichPlot
                case {'raw'}
                    x=expt.X(:,tix);
                    
                case {'norm'}
                    x=expt.Xnorm(:,tix);
                    
                case {'detrend'}
                    x=expt.Xdetrend(:,tix);
                    
                case {'filt'}
                    x=expt.Xfilt(:,tix);
                    doAllPts=true;
                    
                otherwise
                    error([whichPlot, ' is not a supported trace to plot']);
            end
            
            plot(ax,expt.t,x,'k')
            
            if showPts
                for i=1:expt.nS
                    tt=expt.t(expt.segment(i).ix);
                    xx=x(expt.segment(i).ix);
                    pts=expt.points{i}(tix);
                    hold on
                    xPer=interp1(tt,xx,pts.tPer); %this is necessary except for filt
                    plot(ax,pts.tPer,xPer,'bs')
                    if doAllPts
                        plot(ax,pts.tMin,pts.xMin,'r^')
                        plot(ax,pts.tMax,pts.xMax,'rv')
                        plot(ax,pts.tUp,pts.xUp,'bd')
                        plot(ax,pts.tDown,pts.xDown,'bo')
                    end
                    hold off
                end
            end
            
            axis tight
            xlabel('t')
            ylabel('x')
            YLIM=ylim();
            ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
            
            hold on
            for i=2:expt.nS
                plot(ax,expt.segment(i).endpoints(1)*[1,1],ylim(),'g')
            end
            hold off
            
        end
        
        function plot_psd(expt,ax,tix)
            hold off
            for i=1:expt.nS
                plot(ax,expt.f,expt.psd(:,tix,i)); 
%                 plot(ax,expt.f,pow2db(expt.psd(:,tix,i)));
                hold on
            end
            hold off
            % set(gca,'yscale','log')
%             fhi=find(
%             title('One Sided Power Spectral Density');       
            xlabel('frequency')         
            ylabel('power');
            axis tight
            xlim([0,1])
        end
        
        function plot_features(expt,ax,xname,yname,tix)
            
            xname=validatestring(xname,['t','segment',expt.fnames]);
            yname=validatestring(yname,expt.fnames);
            
            
            
            switch xname
                case {'t'}
                    
                case {'segment',''}
                    
                case expt.fnames
                    
                otherwise
                    error(['Invalid name for x-axis: ' xname])
            end
            
        end
        
        function traceKeypress(expt,src,event,whichPlot,showPts)
            tix=src.UserData;
            ax=gca;
            switch(event.Key)
                case {'leftarrow'}
                    if tix>1
                        tix=tix-1;
                        expt.plot_t(ax,whichPlot,tix,showPts)
                    end
                case {'rightarrow'}
                    if tix<expt.nX
                        tix=tix+1;
                        expt.plot_t(ax,whichPlot,tix,showPts)
                    end
            end
            src.UserData=tix;
        end
        
        function psdKeypress(expt,src,event)
            tix=src.UserData;
            ax=gca;
            switch(event.Key)
                case {'leftarrow'}
                    if tix>1
                        tix=tix-1;
                        expt.plot_psd(ax,tix)
                    end
                case {'rightarrow'}
                    if tix<expt.nX
                        tix=tix+1;
                        expt.plot_psd(ax,tix)
                    end
            end
            src.UserData=tix;
        end
    end

end