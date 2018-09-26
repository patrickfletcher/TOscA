%script version of analysis
% interactive plotting at each stage
clear

doPlot=1;

%load data
% filename='/home/pfletcher/Dropbox/CurrentWork/Islets/PKARoscillations/data/#1508_PKAR_11G_18.07.20a.xlsx';
% filename='/home/pfletcher/Dropbox/CurrentWork/Islets/PKARoscillations/data/#1487_PKAR_11G_18.06.22a.xlsx';
filename='/home/pfletcher/Dropbox/CurrentWork/Islets/SatinLabIsletOscillations/data/8G-11G_female RIPCre control_+_+_cre+.xls';
% [filename,path]=uigetfile({'*.xls*'});
% filename=[path,filename];

%load data
[t,XRaw,opts]=loadData(filename);
ix=t<38; t=t(ix); XRaw=XRaw(ix,:);
% ix=t>38; t=t(ix); XRaw=XRaw(ix,:);


%normalize
% XNorm=XRaw;
XNorm=normalizeTraces(t,XRaw,'devmean',[],doPlot);
% XNorm=normalizeTraces(t,XRaw,'devtrend',{'gaussian',15},doPlot);

%detrend
% XDT=XNorm;
XDT=detrendTraces(t,XNorm,'gaussian',20,doPlot);

XDT=XDT-mean(XDT,1);

%filter
XFilt=filterTraces(t,XDT,'gaussian',1,doPlot);

%% plot psd
fs=1/mode(diff(t));
% [PN,fN]=powerSpectrum(XNorm,fs);
% [PDT,fDT]=powerSpectrum(XDT,fs);
% [PF,fF]=powerSpectrum(XFilt,fs);
% [PN,fN]=pwelch(XNorm,[],[],4096,fs);
% [PDT,fDT]=pwelch(XDT,[],[],4096,fs);
% [PF,fF]=pwelch(XFilt,[],[],4096,fs);
wind=window(@hamming,length(t));
[PN,fN]=periodogram(XNorm,wind,4096,fs);
[PDT,fDT]=periodogram(XDT,wind,4096,fs);
[PF,fF]=periodogram(XFilt,wind,4096,fs);
%%
tix=9;
figure(1)
plot(fN,pow2db(PN(:,tix)),fDT,pow2db(PDT(:,tix)),fF,pow2db(PF(:,tix)));
%%

%final rescaling

%feature detection
