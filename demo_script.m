%script version of analysis
% interactive plotting at each stage
clear

doPlot=1;

%load data
filename='/home/pfletcher/Dropbox/CurrentWork/Islets/PKARoscillations/data/#1508_PKAR_11G_18.07.20a.xlsx';
% filename='/home/pfletcher/Dropbox/CurrentWork/Islets/PKARoscillations/data/#1487_PKAR_11G_18.06.22a.xlsx';
% filename='/home/pfletcher/Dropbox/CurrentWork/Islets/SatinLabIsletOscillations/data/8G-11G_female RIPCre control_+_+_cre+.xls';
% [filename,path]=uigetfile({'*.xls*'});
% filename=[path,filename];

%load data
[t,XRaw,opts]=loadData(filename);

%normalize
XNorm=normalizeTraces(t,XRaw,'devmean',[],doPlot);
% XNorm=normalizeTraces(t,XRaw,'devtrend',{'gaussian',15},doPlot);

%detrend
% XDT=XNorm;
XDT=detrendTraces(t,XNorm,'gaussian',15,doPlot);

%filter
XFilt=filterTraces(t,XDT,'gaussian',1,doPlot);

%final rescaling

%feature detection
