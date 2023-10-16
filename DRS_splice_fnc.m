function [full_wv, full_DRS] = DRS_splice_fnc(DRS_VIS, DRS_NIR, VIS_wv, NIR_wv, bkgd, ref, tissue_type, NumOfMeas, save_flags)
% Combine and calibrate DRS measurements by VIS and NIR spectrometers from OceanOptics.
% Celina L. Li, June 2021.

% Inputs:
% DRS_VIS = raw DRS - VIS
% DRS_NIR = raw DRS - NIR
% VIS_wv = VIS wavelength lambda in nm
% NIR_wv = NIR wavelength lambda in nm
% bkgd - background measurements in {VIS NIR}.
% ref - reference measurements in {VIS NIR}.
% tissue_type - the tissue type in string e.g. 'boneMarrow' .
% NumOfMeas - number of locations measured.
% save_flag = [1,1,1] indicates TO SHOW PLOT, SAVE FIGURE, SAVE PROC DATA.
% save_flag = [0,0,0] indicates TO NOT SHOW PLOT, NOT SAVE FIGURE, NOT SAVE PROC DATA.

% Linear interpolation of VIS & NIR spectral measurements was adapted from
% Nogueira, M. S. et al (2021). DOI: 10.1038/s41598-020-79517-2

% ===========================================================================
% Calibration - Background & Reference

NIR_avg_cal = []; VIS_avg_cal = [];

for ii = 1:NumOfMeas % number of locations measured
    
    % NIR
    % calculate calibrated average; data format 512 x #of locations
    NIR_avg_cal_holder = ( DRS_NIR(:,ii)-bkgd{2} )./( ref{2}-bkgd{2} );
    NIR_avg_cal = cat( 2, NIR_avg_cal, NIR_avg_cal_holder );

    % VIS
    % calculate calibrated average; data format 512 x #of locations
    VIS_avg_cal_holder = ( DRS_VIS(:,ii)-bkgd{1} )./( ref{1}-bkgd{1} );
    VIS_avg_cal = cat( 2, VIS_avg_cal, VIS_avg_cal_holder );

end

%Combine VIS & NIR

VIS_ind_init = find(VIS_wv > 355);
VIS_ind = find(1095 < VIS_wv & VIS_wv < 1120);
NIR_ind = find(1095 < NIR_wv & NIR_wv < 1120);

VIS_wv_new = VIS_wv( VIS_ind_init(1):VIS_ind(1)-1 );
VIS_avg_cal_new = VIS_avg_cal(VIS_ind_init(1):VIS_ind(1)-1,:);

VIS_wv_com = VIS_wv( VIS_ind(1):VIS_ind(end) );
VIS_avg_cal_com = VIS_avg_cal( VIS_ind(1):VIS_ind(end),: );

NIR_wv_new = NIR_wv( NIR_ind(end)+1:end);
NIR_avg_cal_new = NIR_avg_cal( NIR_ind(end)+1:end,: );

NIR_wv_com = NIR_wv( NIR_ind(1):NIR_ind(end) );
NIR_avg_cal_com = NIR_avg_cal( NIR_ind(1):NIR_ind(end),: );

% Spectral Overlapping Interpolation

n = 100 + 1; 

NIR_spec_wt = linspace(0, 1, n);
VIS_spec_wt = (fliplr(NIR_spec_wt))';

overlap_wv = ( linspace( VIS_wv(VIS_ind(1)),NIR_wv(NIR_ind(end)),n ) )';

full_wv = [VIS_wv_new; overlap_wv; NIR_wv_new];

for i = 1:NumOfMeas
    VIS_spec_interp(:,i) = interp1( VIS_wv_com,VIS_avg_cal_com(:,i),overlap_wv,'spline' );
    NIR_spec_interp(:,i) = interp1( NIR_wv_com,NIR_avg_cal_com(:,i),overlap_wv,'spline' ); 
   
    common_spec(:,i) = VIS_spec_wt.*VIS_spec_interp(:,i)+( sum(VIS_spec_interp(:,i))/sum(NIR_spec_interp(:,i)) )*NIR_spec_wt'.*NIR_spec_interp(:,i);
    
    full_DRS(:,i) = [VIS_avg_cal_new(:,i);common_spec(:,i);(sum(VIS_spec_interp(:,i))/sum(NIR_spec_interp(:,i)))*NIR_avg_cal_new(:,i)];
end

% Plot

if save_flags(1) == 1
    fig = figure;
    plot(full_wv, full_DRS/max(full_DRS,[],'all'));

    title([tissue_type ' DRS'],'FontSize',14,'FontWeight','bold');
    xlabel('Wavelength (nm)','FontSize',12,'FontWeight','bold');
    ylabel('Normalized Intensity (a.u.)','FontSize',12,'FontWeight','bold');
    xlim([300 2000]);
    ylim([0 1]);

    x0=10;
    y0=10;
    width=850;
    height=500;
    set(gcf,'position',[x0,y0,width,height]);
    set(gca,'Xtick',300:100:2000)
    grid on   
end
    
% Save figure
if save_flags(2) == 1
    
    sprintf('Please select save path for figure')
    save_dir_fig = uigetdir;
    
    fig_name_1 = strcat(tissue_type, '.jpg');
    %fig_name_2 = strcat(tissue_type, '.fig');
    
    if isfile( fullfile(save_dir_fig, fig_name_1) )
    else
    saveas( fig, fullfile(save_dir_fig, fig_name_1) );
    end
    
    %if isfile( fullfile(save_dir_fig, fig_name_2) )
    %else
    %saveas( fig, fullfile(save_dir_fig, fig_name_2) );
    %end
    
end

% Save proc. data

full_wv = full_wv.';
full_DRS = full_DRS.';

if save_flags(3) == 1
    
    sprintf('Please select save path for processed data')
    save_dir_data = uigetdir;
    
    file_name = strcat(tissue_type, '_DRS.mat');
    
    if isfile( fullfile(save_dir_data, file_name) )
    else
        save(fullfile(save_dir_data, file_name),'full_DRS');
    end
    
    if isfile( fullfile(save_dir_data, '\', tissue_type, '_DRS.csv') )
    else    
        writematrix(full_DRS,strcat(save_dir_data, '\', tissue_type, '_DRS.csv'));
    end
    
    if isfile( fullfile(save_dir_data, 'full_DRS_wv.mat') )
    else
        save(fullfile(save_dir_data, 'full_DRS_wv.mat'), 'full_wv');
    end
    
    if isfile( fullfile(save_dir_data, '\full_DRS_wv.csv') )
    else
        writematrix(full_wv,strcat(save_dir_data, '\full_DRS_wv.csv'));
    end
    
end