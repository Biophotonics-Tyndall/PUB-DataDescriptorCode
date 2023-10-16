%% Calibration and Pre-Processing of raw DRS measurements
% 1. Calibration using (S_raw - S_bkgd)/(S_ref - S_bkgd); using
% DRS_splice_fnc
% 2. Splice calibrated VIS & NIR DRS together; using
% DRS_splice_fnc

% Correct a systematic 'defect' between 1250 - 1300 nm using moving median

clear all 
close all
clc

dir_path = uigetdir; %Select Main Data Folder

cd(dir_path)  %change to that directory

%[p,dir_name,c]=fileparts(dir_path);

%% Load files

cd(dir_path);
cd 'raw_mat';

% Background
load('NIR_bkgd_avg.mat');
load('VIS_bkgd_avg.mat');

% Reference
load('NIR_ref_avg.mat');
load('VIS_ref_avg.mat');

% Wavelengths
load('NIR_wv.mat');
load('VIS_wv.mat');

%% Plot Background

% NIR
figure('Name','NIR Bkgd','NumberTitle','off');
hold on
plot( NIR_wv, NIR_bkgd_avg(2,:) );
hold off
title('NIR Bkgd','FontSize',14,'FontWeight','bold');
xlabel('Wavelength (nm)','FontSize',12,'FontWeight','bold');
ylabel('Intensity (a.u.)','FontSize',12,'FontWeight','bold');

x0=10;
y0=10;
width=850;
height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'Xtick',300:100:2000)
grid on

% VIS
figure('Name','VIS Bkgd','NumberTitle','off');
hold on
plot( VIS_wv, VIS_bkgd_avg(2,:) );
hold off
title('VIS Bkgd','FontSize',14,'FontWeight','bold');
xlabel('Wavelength (nm)','FontSize',12,'FontWeight','bold');
ylabel('Intensity (a.u.)','FontSize',12,'FontWeight','bold');

x0=10;
y0=10;
width=850;
height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'Xtick',300:100:2000)
grid on

%% Plot Reference (Reflectance Standard)

% NIR
figure('Name','NIR Ref','NumberTitle','off');
hold on
plot( NIR_wv, NIR_ref_avg(2,:) );
hold off
title('NIR Ref','FontSize',14,'FontWeight','bold');
xlabel('Wavelength (nm)','FontSize',12,'FontWeight','bold');
ylabel('Intensity (a.u.)','FontSize',12,'FontWeight','bold');

x0=10;
y0=10;
width=850;
height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'Xtick',300:100:2000)
grid on

% VIS
figure('Name','VIS ref','NumberTitle','off');
hold on
plot( VIS_wv, VIS_ref_avg(2,:) );
hold off
title('VIS Ref','FontSize',14,'FontWeight','bold');
xlabel('Wavelength (nm)','FontSize',12,'FontWeight','bold');
ylabel('Intensity (a.u.)','FontSize',12,'FontWeight','bold');

x0=10;
y0=10;
width=850;
height=500;
set(gcf,'position',[x0,y0,width,height]);
set(gca,'Xtick',300:100:2000)
grid on

%% Example - Cartilage

cd(dir_path);
cd 'raw_mat';

% Load Raw DRS Data - VIS & NIR
load('cartilage_DRS_NIR.mat');
load('cartilage_DRS_VIS.mat');

tissue_type = 'cartilage';
bkgd = {VIS_bkgd_avg(2,:)', NIR_bkgd_avg(2,:)'};
ref = {VIS_ref_avg(2,:)', NIR_ref_avg(2,:)'};
NumOfMeas = length(NIR_avg(:,1));
save_flags = [1,1,1];

% To splice VIS & NIR
[full_wv, full_DRS] = DRS_splice_fnc(VIS_avg', NIR_avg', VIS_wv', NIR_wv', bkgd, ref, tissue_type, NumOfMeas, save_flags);

%% Spectral Correction - Smoothing using Moving Median

idx_wv = find(full_wv<=1850);
all_spec = full_DRS(:,idx_wv);
full_wv = full_wv(:,idx_wv);

smooth_idx = find(full_wv>=1250&full_wv<=1300);
smooth_spec = smoothdata(all_spec(:,smooth_idx),2,'movmedian',2.2);
all_spec(:,smooth_idx) = smooth_spec;

