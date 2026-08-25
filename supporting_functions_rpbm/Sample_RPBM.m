clear all; close all; clc;
%% Sample RPBM
% Sample fitting and diffusion parameter estimation from the random 
% permeable barrier model (RPBM) on D(t). This sample script calls
% get_Dt_RPBM to estimate parameters tau and zeta. Afterwards, RPBM_process
% is used to parse these parameters and output relevant biophysical. 
%
% Novikov   et al. Nature Physics 7:508 (2011)
% Fiermeans et al. NMR in Biomed 30:e3612 (2017)
%
% (c) Gregory Lemberskiy 2017 gregorylemberskiy@gmail.com

%% Input - Sample Soleus Muscle D(t) following DTI Tensor Fitting
% Soleus ROI from a DWI acquisition using 3x3x5 mm^3 voxels, 20 directions
% at b=500 s/mm^2 and 3 b=0. Outlier rejection was performed using IRESTORE
% and tensor estimation was performed using weighted linear least squares.
%
% Chang   et al. IRESTORE. Magn Reson Med 68: 1654?1663 (2012) 
% Veraart et al. WLLS    .  Neuroimage     81:  335-46   (2013)

time = [40,100,200,500,800];                              %Diffusion Time
DL   = [1.62518,1.62594,1.44837,1.55005,1.450938];   %Longitudinal Diffusivity
DR   = [1.02629,0.93104,0.77363,0.71837,0.554827];   %Radial Diffusivity

uDL=[0.25175,0.18469,0.21147,0.22214,0.18235];
uDR=[0.18904,0.14447,0.11496,0.14323,0.11731];

UseWeights=1;

%% RPBM FITTING

FIX_D0=0;
if UseWeights==1; W=1./(uDR.^2); else; W=1; end     %Weights from std of Dr(t)
% W = [1 1 1 1 1];

%% Fitting by Fixing D0 to L1(t~inf)
% In principle, D0 should be equal to L1 at long times. However, in the
% presense of noisy data this the eigenvalues will repel from one another
% causing a systematic bias on parameter estimation.
%
% Conversely if the structure of interest is highly restricted (small
% diameters), diffusion times may be too long to accurately estimate D0.
Dfix = mean(DL(time>200));                                                 %Fixing D0 to L1(t~inf)
uDfix = std(DL(time>200));  

opt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[Inf,Inf],...
    'Startpoint',[300,2],...
    'Robust','bisquare',...
    'Weights',W);
fitDapp = fittype('D*get_Dt_RPBM(x/tau,zeta)','problem','D','options',opt);
[fitresult_F,gof_F] = fit(time',DR',fitDapp,'problem',double(Dfix));
[RPBM_F,uRPBM_F]=RPBM_Process([Dfix,coeffvalues(fitresult_F)],[uDfix,uDfix;confint(fitresult_F)']);%Calculate Various RPBM Parameters and save in structure
txn=linspace(1,2000,2000);
RPBM_FIX=squeeze(fitDapp(fitresult_F.tau,fitresult_F.zeta,fitresult_F.D,txn));

%% Fitting by Varying D0
% If there is sufficient SNR and sampling of t, fitting over 3 parameters
% would be the ideal approach. Given your acquisition, I recommend trying
% both approaches [fitting/fixing D0] to see if they are in agreement.

Dfix = mean(DL(time>200));                                                 %Starting Value at D0~L1(t~inf)
opt = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'Upper',[Inf,Inf,Inf],...
    'Startpoint',[Dfix,300,2],...
    'Robust','bisquare',...
    'Weights',W);
fitDapp = fittype('D*get_Dt_RPBM(x/tau,zeta)','options',opt);
[fitresult_V,gof] = fit(time',DR',fitDapp);
[RPBM_V,uRPBM_V]=RPBM_Process([coeffvalues(fitresult_V)],[confint(fitresult_V)']);%Calculate Various RPBM Parameters and save in structure

txn=linspace(1,2000,2000);
RPBM_VARY=squeeze(fitDapp(fitresult_V.D,fitresult_V.tau,fitresult_V.zeta,txn));

%% Plotting Fits
% Much of the transient time-dependence of D(t) happens at the shortest
% diffusion times. It is often valuable to inspect D(t) on a semilog scale.

figure(1)
subplot(1,2,1)
plot(time,DL,'rs','markerfacecolor','r','markersize',10); hold on;
plot(time,DR,'ro','markerfacecolor','r','markersize',10);
x1(1)=plot(txn,RPBM_FIX,'r--','linewidth',2);
x1(2)=plot(txn,RPBM_VARY,'b--','linewidth',2);
xlabel('Time [ms]','fontsize',20)
ylabel('D_{||}, D_\perp','fontsize',20)
set(gca,'linewidth',2,'fontweight','bold')
pbaspect([1 1 1])
ylim([0.8 1.9])
xlim([0 800])

subplot(1,2,2)
plot(time,DL,'rs','markerfacecolor','r','markersize',10); hold on;
plot(time,DR,'ro','markerfacecolor','r','markersize',10);
plot(txn,RPBM_FIX,'r--','linewidth',2)
plot(txn,RPBM_VARY,'b--','linewidth',2)
xlabel('Time [ms]','fontsize',20)
ylabel('D_{||}, D_\perp','fontsize',20)
set(gca,'linewidth',2,'fontweight','bold')
set(gca,'xscale','log')
set(gca,'xtick',[1 10 100 500 2500 10000])
title('Semilog X','fontsize',20)
pbaspect([1 1 1])
ylim([0.8 1.9])
xlim([txn(1) txn(end)])

