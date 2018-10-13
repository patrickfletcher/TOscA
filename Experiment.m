classdef Experiment < handle
    %class to represent a single experiment containing >=1 timeseries that share the same time sample points. For
    %example, a fluorescence videomicroscopy experiment with multiple ROIs in a field of view. 

    properties
        name=''
        date=''
        sex=''
        animalCondition='' %eg. WT vs KO
        
        appNames={''} %eg. '5G', 'Dz',..
        appTimes=[] %list of start times for each condition
        tIntervals=[]
        
        filename=''
        notes={}
        
        t
        dt
        fs
        
        X
        Xnorm
        Xdetrend
        Xfilt
        
        include %set element to zero to exclude a trace
        nX
        
        group
        Xg %average traces with same group (always use Rdetrend?)
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
            
            %instead to interpolate on a fixed grid with given DT
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
            
            % time intervals
            if ~exist('appTimes','var')
                appTimes=[];
            end
            appTimes=[expt.t(1),appTimes(:)',expt.t(end)];
            for i=1:length(appTimes)-1
                expt.tIntervals(i,:)=[appTimes(i),appTimes(i+1)];
            end
            expt.appTimes=appTimes(1:end-1);
            
            % omit fixed amount of time at beginning/end of intervals?
%             tOmit=0;
%             expt.nOmit=round(expt.tOmit/expt.dt); 

            expt.include=true(1,expt.nX);
        end
        
        %save results. If required metadata is not set, force user to enter it now
        function save()
        end
        
        %simple wrappers for preprocessing functions for now
        %TODO: forced order? if so, simply do a preprocess?
        function normalize(expt,method,methodPar,doPlot)
            expt.Xnorm=normalizeTraces(expt.t,expt.X,method,methodPar,doPlot);
            expt.normMethod=method;
            expt.normParam=methodPar;
        end
        
        function detrend(expt,method,methodPar,doPlot)
            if isempty(expt.Xnorm)
                expt.Xnorm=expt.X;
            end
            expt.Xdetrend=detrendTraces(expt.t,expt.Xnorm,method,methodPar,doPlot);
            expt.trendMethod=method;
            expt.trendParam=methodPar;
        end
        
        function filter(expt,method,methodPar,doPlot)
            if isempty(expt.Xnorm)
                expt.Xnorm=expt.X;
            end
            if isempty(expt.Xdetrend)
                expt.Xdetrend=expt.X;
            end
            expt.Xfilt=filterTraces(expt.t,expt.Xdetrend,method,methodPar,doPlot);
            expt.filterMethod=method;
            expt.filterParam=methodPar;
        end
        
    end

end