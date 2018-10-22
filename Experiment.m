classdef Experiment < handle
    %class to represent a single experiment containing >=1 timeseries that share the same time sample points. For
    %example, a fluorescence videomicroscopy experiment with multiple ROIs in a field of view. 

    properties
        name=''
        date=''
        sex=''
        condition='' %eg. WT vs KO
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
        
        thrFrac=[0.55,0.45]
        delta=0.25
        Tbig=15
        Tsmall=4
        
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
        
        function setGroup(expt,groupvar)
            expt.group=groupvar;
            expt.nG=unique(groupvar);
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
        
        
        %Analysis functions - per interval
        function compute_features(expt,method,methodPar)
            expt.detect_periods(method,methodPar);
            expt.periodogram();
            
            expt.fnames=fieldnames(expt.features{1})';
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
            if ~isempty(expt.psd)
                       
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
        end
        
        function plotFeatures(expt,ax,xname,yname,tix)
            %xname={t,segment,fnames}
            %yname={fnames}
                
            if ~isempty(expt.features)
                
            if isempty(ax)
                ax=gca;
            end
            
            if isempty(xname)
                xname='segment';
            end
            
            xname=validatestring(xname,['t','segment',expt.fnames]);
            yname=validatestring(yname,expt.fnames);
            
            if ~exist('tix','var')||isempty(tix)
                tix=1;
                figID=gcf; %newfig?
                figID.KeyPressFcn={@expt.featureKeypress,xname,yname};
                figID.UserData=tix; %store the current trace ID with figure
                ax=gca;
            end
            
            expt.plot_features(ax,xname,yname,tix);
            
            end
        end
        
    end
    
    methods (Access=private)
        
        
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
                
                %also get xvalues for all X types
                %independently max/min per period for non-filt traces (detrend at least)
                
                %get average features per segment per trace
                Tmean=arrayfun(@(x)mean(x.period),feats); Tmean=num2cell(Tmean);
                Tstd=arrayfun(@(x)std(x.period),feats); Tstd=num2cell(Tstd);
                APDmean=arrayfun(@(x)mean(x.APD),feats); APDmean=num2cell(APDmean);
                APDstd=arrayfun(@(x)std(x.APD),feats); APDstd=num2cell(APDstd);
                PFmean=arrayfun(@(x)mean(x.PF),feats); PFmean=num2cell(PFmean);
                PFstd=arrayfun(@(x)std(x.PF),feats); PFstd=num2cell(PFstd);
                Amean=arrayfun(@(x)mean(x.amp),feats); Amean=num2cell(Amean);
                Astd=arrayfun(@(x)std(x.amp),feats); Astd=num2cell(Astd);
                
                [feats.Tmean]=Tmean{:};
                [feats.Tstd]=Tstd{:};
                [feats.APDmean]=APDmean{:};
                [feats.APDstd]=APDstd{:};
                [feats.PFmean]=PFmean{:};
                [feats.PFstd]=PFstd{:};
                [feats.Amean]=Amean{:};
                [feats.Astd]=Astd{:};
                
                expt.points{i}=pts;
                expt.features{i}=feats;
            end
%             expt.fnames=fieldnames(feats)';
        end
        
        function periodogram(expt)
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                [PSD,F,Pmax,fmax]=powerSpectrum(expt.Xfilt(thisIx,:),expt.fs);
                Tpsd=1./fmax;
                expt.psd(:,:,i)=PSD;
                
                Tpsd=num2cell(Tpsd);
                fmax=num2cell(fmax);
                Pmax=num2cell(Pmax);
                [expt.features{i}.Tpsd]=Tpsd{:};
                [expt.features{i}.fmax]=fmax{:};
                [expt.features{i}.Pmax]=Pmax{:};
            end
            expt.f=F;
%             expt.fnames=fieldnames(expt.features{i})'; 
        end
        
        function plot_t(expt,ax,whichPlot,tix,showPts)
            
            axes(ax)
            cla
            hold(ax,'on');
            
            doAllPts=false;
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
                    doAllPts=true;
                    
                otherwise
                    error([whichPlot, ' is not a supported trace to plot']);
            end
            
            plot(ax,expt.t,x,'k')
            if ~isempty(x2)
            plot(ax,expt.t,x2)
            end
            
            if showPts
                for i=1:expt.nS
                    tt=expt.t(expt.segment(i).ix);
                    xx=x(expt.segment(i).ix);
                    pts=expt.points{i}(tix);
                    
                    if length(pts.tPer)>1
                        
                    xPer=interp1(tt,xx,pts.tPer); %this is necessary except for filt
                    plot(ax,pts.tPer,xPer,'bs')
                    if doAllPts
                        plot(ax,pts.tMin,pts.xMin,'r^')
                        plot(ax,pts.tMax,pts.xMax,'rv')
                        plot(ax,pts.tUp,pts.xUp,'bd')
                        plot(ax,pts.tDown,pts.xDown,'bo')
                        
                        tupdwn=[pts.tUp(1:length(pts.tDown)); pts.tDown];
                        ythresh=[1;1]*expt.features{i}(tix).pthresh;
                        plot(ax,tupdwn,ythresh,'b-')
                    end
                    end
                end
            end
            
            axis tight
            xlabel('t')
            ylabel(['X',whichPlot])
            YLIM=ylim(ax);
            ylim(ax,[YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
            
            for i=2:expt.nS
                plot(ax,expt.segment(i).endpoints(1)*[1,1],ylim(),'g')
            end
            hold(ax,'off')
            
        end
        
        function plot_psd(expt,ax,tix)
            
            axes(ax)
            cla
            hold(ax,'on');
            
            for i=1:expt.nS
                plot(ax,expt.f,expt.psd(:,tix,i)); 
%                 plot(ax,expt.f,pow2db(expt.psd(:,tix,i)));
            end
            % set(gca,'yscale','log')
%             fhi=find(
%             title('One Sided Power Spectral Density');       
            xlabel('frequency')         
            ylabel('power');
            axis tight
            xlim(ax,[0,0.5])
            if expt.nS>1
                legend({expt.segment.name})
            end
            hold(ax,'off')
        end
        
        function plot_features(expt,ax,xname,yname,tix)
            
            axes(ax)
            cla
            hold(ax,'on');
            
            
            %x,y - highlight point for scalar features
            %xx,yy - distribution of features (either across all traces for scalar features, or distribution within a
            %trace for feature distributions)
            
            xJitter=0.05;
           
            HX=[];HY=[];
            XX={};YY={};
            
            
            scalarFeatures=false;
            colorBySegment=false;
            for i=1:expt.nS
                        
                %TODO: more robust way? what if that particular trace had 1 or 0 periods detected
                x=[];
                y=expt.features{i}(tix).(yname);
                if isempty(y)
                    return
                end
                
                switch xname
                    case {'t'}
                        xx=expt.points{i}(tix).tMax;
                        if isscalar(y) %must be a distribution feature
                            error('Only per-trace feature distributions may be plotted vs time')
                        else
                            yy=y;
                        end
                        
                        
                    case {'segment'}
                        if isscalar(y) %plot all trace data as line seg, highlight point tix
%                             xx=i+xJitter*randn(1,expt.nX);
                            xx=i*ones(1,expt.nX);
                            yy=[expt.features{i}.(yname)];
                            scalarFeatures=true;
                            HX(end+1,1)=i;
                            HY=[HY;y(:)];
                        else
                            xx=i+xJitter*randn(1,length(y));
                            yy=y;
                        end
                        

                    case expt.fnames
                        x=expt.features{i}(tix).(xname);
                        
                        if isscalar(x) %plot all trace data, highlight point tix
                            if isscalar(y) %plot all trace data, highlight point tix
                                xx=[expt.features{i}.(xname)];
                                yy=[expt.features{i}.(yname)];
                                scalarFeatures=true;
                                HX=[HX;x(:)];
                                HY=[HY;y(:)];
                            else
                                error('x and y feature type must match')
                            end
                        else
                            if isscalar(y) %plot all trace data, highlight point tix
                                error('x and y feature type must match')
                            else
                                xx=x;
                                yy=y;
                                colorBySegment=true;
                            end
                        end

                    otherwise
                        error(['Invalid name for x-axis: ' xname])
                end
                
                XX=[XX;{xx}];
                YY=[YY;{yy}];

            end
                            
            if scalarFeatures
                
                XX=cell2mat(XX);
                YY=cell2mat(YY);
                
                if expt.nS>1
                    plot(ax,XX,YY,'k-')
                    plot(ax,HX,HY,'k-','linewidth',1.5)
                end
                h2=plot(ax,XX',YY','o');
                for i=1:expt.nS
                    h4=plot(ax,HX(i),HY(i),'o');
                    h4.Color=h2(i).Color;
                    h4.MarkerFaceColor=h2(i).Color;
                end
                if expt.nS>1
                    legend(h2,{expt.segment.name})
                end
            else
                for i=1:expt.nS
                    hl(i)=plot(ax,XX{i},YY{i},'o');
                end
                if colorBySegment
                    legend(hl,{expt.segment.name})
                else
                    for i=1:expt.nS
                        hl(i).Color='k';
                    end
                end
            end
            
            axis tight
            YLIM=ylim(ax);
            ylim(ax,[YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
            
            switch xname
                    case {'t'}
                        xlabel('t')
                        for i=1:expt.nS
                            hl(i).LineStyle='-';
                        end
                        
                        for i=2:expt.nS
                            plot(ax,expt.segment(i).endpoints(1)*[1,1],ylim(),'g')
                        end
                        
                    case {'segment'}
                        xticks(ax,1:expt.nS)
                        xticklabels(ax,{expt.segment.name})
                        xlim(ax,[0.5,expt.nS+0.5])
                        xlabel(ax,'segment')
                        
                    case expt.fnames
                        xlabel(ax,xname)
            end
            
            ylabel(ax,yname)
            
            hold(ax,'off')
            
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
        
        function featureKeypress(expt,src,event,xname,yname)
            tix=src.UserData;
            ax=gca;
            switch(event.Key)
                case {'leftarrow'}
                    if tix>1
                        tix=tix-1;
                        expt.plot_features(ax,xname,yname,tix);
                    end
                case {'rightarrow'}
                    if tix<expt.nX
                        tix=tix+1;
                        expt.plot_features(ax,xname,yname,tix);
                    end
            end
            src.UserData=tix;
        end
    end

end