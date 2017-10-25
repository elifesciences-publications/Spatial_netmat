clear all; close all; clc

% % Set up paths
% Dir = '/home/fs0/janineb/scratch/';
% path(path,fullfile(Dir,'matlab'))
% path(path,fullfile(Dir,'HCP','CCA','files'))
% path(path,fullfile(Dir,'matlab','cifti-matlab'))
% path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
% addpath ~steve/NETWORKS/FSLNets;
% path(path,fullfile(Dir,'HCP','netmat_simulations','palm-alpha106'));
% 
% % Load input data
% load('Results/input_CCA_table1.mat','Glasser_amp','Glasser_maps');
% load('Results/input_CCA_Glasser.mat','FccaORIG_imputed','PccaORIG_imputed');
% load('Results/input_CCA_table1_PFM_441.mat');
% load('Results/input_CCA_table1_fractional_area.mat');
% 
% % Set number of permutations
Nperm = 100000;
% 
% % Create output variables
% CCAall = zeros(9,5);
% NULLall = zeros(Nperm,9);
% Uall = zeros(441,9); Vall = zeros(441,9);
% 
% % Run all CCAs
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(Glasser_amp,Nperm);
% CCAall(1,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,1) = Null(:,1); Uall(:,1) = U(:,1); Vall(:,1) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(Glasser_maps,Nperm);
% CCAall(2,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,2) = Null(:,1); Uall(:,2) = U(:,1); Vall(:,2) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(FccaORIG_imputed,Nperm);
% CCAall(3,:) = [U(1) V(1) R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,3) = Null(:,1); Uall(:,3) = U(:,1); Vall(:,3) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(PccaORIG_imputed,Nperm);
% CCAall(4,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,4) = Null(:,1); Uall(:,4) = U(:,1); Vall(:,4) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(PFM441_amplitude,Nperm);
% CCAall(5,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,5) = Null(:,1); Uall(:,5) = U(:,1); Vall(:,5) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(PFM441_spatial,Nperm);
% CCAall(6,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,6) = Null(:,1); Uall(:,6) = U(:,1); Vall(:,6) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(PFM441_Fnetmat,Nperm);
% CCAall(7,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,7) = Null(:,1); Uall(:,7) = U(:,1); Vall(:,7) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(PFM441_Pnetmat,Nperm);
% CCAall(8,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,8) = Null(:,1); Uall(:,8) = U(:,1); Vall(:,8) = V(:,1);
% 
% [~,U,V,R,Ruuica,Iuuica,P,Puuica,Null] = S_CCA(Fractional_area,Nperm);
% CCAall(9,:) = [R(1) Ruuica(1) Iuuica(1) P(1) Puuica(1)];
% NULLall(:,9) = Null(:,1); Uall(:,9) = U(:,1); Vall(:,9) = V(:,1);
% 
% % Save outputs
% save('Results/CCA_MMP_PFM.mat','NULLall','CCAall','Uall','Vall','-v7.3');

% Look at significance & significant difference across all:
Names = {'MMP amplitude','MMP spatial','MMP full netmat','MMP partial netmat',...
    'PFM amplitude','PFM spatial','PFM full netmat','PFM partial netmat',...
    'Fractional area'};
load('Results/CCA_MMP_PFM.mat');
Pthres = prctile(NULLall,95); figure; plot(Pthres);
DIFFthres = zeros(9*8/2,1); I = 1;
DIFFsig = zeros(9,9);
for x = 1:9
    for y = 1:9
        if y>x
            null_temp = NULLall(:,x)-NULLall(:,y);
            D = CCAall(x,1)-CCAall(y,1);
            DIFFthres(I) = prctile(null_temp,95);
            DIFFsig(x,y) = (1+sum(null_temp(2:end,1)>=D))/Nperm; DIFFsig(y,x) = DIFFsig(x,y);
            I = I+1;
        end
    end
end
clear x y I null_temp
figure; plot(DIFFthres)

% Correlate all Us
DIFFsig(eye(9)==1) = 1;
Ucorr = corr(Uall);
figure; imagesc(abs(Ucorr),[0 1]); 
hold on; [x,y] = ind2sub([9 9],find(DIFFsig<0.05)); text(x-0.1,y+0.1,'*','color','k','fontsize',20);
set(gca,'xtick',1:9,'xticklabel',Names,'ytick',1:9,'yticklabel',Names)
set(gcf,'Position',[100 100 1100 450],'PaperPositionMode','auto')
