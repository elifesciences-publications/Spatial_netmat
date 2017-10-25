clear all; close all; clc

% Paths
addpath ~steve/NETWORKS/FSLNets;
DRdir = '/vols/Scratch/janineb/HCP/netmat_simulations';
Ncomponents = 15;
Nkeep = 100;
Nperm = 0; 
addpath /vols/Scratch/janineb/HCP/CCA/files/
load('/vols/Scratch/janineb/HCP/CCA/files/Permutation_task_100000.mat','conf')
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/

% Netmats:
fprintf('Loading ts\n');
ts_grp = nets_load(fullfile(DRdir,'Results','DRtask_groupmaps'),0.72,1,1);
Pgrp = nets_netmats(ts_grp,-1,'ridgep',0.01);
Fgrp = nets_netmats(ts_grp,-1,'corr');
ts_sub = nets_load(fullfile(DRdir,'Results','DRtask_subjectmaps'),0.72,1,1);
Psub = nets_netmats(ts_sub,-1,'ridgep',0.01);
Fsub = nets_netmats(ts_sub,-1,'corr');
ts_subM = nets_load(fullfile(DRdir,'Results','DRtask_subjectmaps','masked_ts'),0.72,1,1);
PsubM = nets_netmats(ts_subM,-1,'ridgep',0.01);
FsubM = nets_netmats(ts_subM,-1,'corr');

% Paired comparison
Bonf = 0.05/(Ncomponents*(Ncomponents-1)/2);
[hF, pF, ~, statsF] = ttest(Fgrp,Fsub);
[hP, pP, ~, statsP] = ttest(Pgrp,Psub);
fprintf('Number of significant differences in full = %d\n',length(find(pF<Bonf)));
fprintf('Number of significant differences in partial = %d\n',length(find(pP<Bonf)));

% CCA
A = Fgrp; A(166,:) = [];
A = demean(A); A = demean(A-conf*(pinv(conf)*A)); [A,~,~]=nets_svds(A,Nkeep);
CCAinput.Fgrp = A;
A = Fsub; A(166,:) = []; 
A = demean(A); A = demean(A-conf*(pinv(conf)*A)); [A,~,~]=nets_svds(A,Nkeep); 
CCAinput.Fsub = A;
A = Pgrp; A(166,:) = []; 
A = demean(A); A = demean(A-conf*(pinv(conf)*A)); [A,~,~]=nets_svds(A,Nkeep); 
CCAinput.Pgrp = A;
A = Psub; A(166,:) = []; 
A = demean(A); A = demean(A-conf*(pinv(conf)*A)); [A,~,~]=nets_svds(A,Nkeep); 
CCAinput.Psub = A;
A = FsubM; A(166,:) = []; 
A = demean(A); A = demean(A-conf*(pinv(conf)*A)); [A,~,~]=nets_svds(A,Nkeep); 
CCAinput.FsubM = A;
A = PsubM; A(166,:) = []; 
A = demean(A); A = demean(A-conf*(pinv(conf)*A)); [A,~,~]=nets_svds(A,Nkeep); 
CCAinput.PsubM = A;

% Group means
[Z_Fgrp,M_Fgrp]=nets_groupmean(Fgrp,0);
[Z_Fsub,M_Fsub]=nets_groupmean(Fsub,0);
[Z_Pgrp,M_Pgrp]=nets_groupmean(Pgrp,0);
[Z_Psub,M_Psub]=nets_groupmean(Psub,0);
[dpRSN,yyRSN] = nets_hierarchy(M_Fgrp,M_Fgrp,1:Ncomponents,'');
figure; 
subplot(2,2,1); imagesc(M_Fgrp(dpRSN,dpRSN),[-1 1]); title('Full netmat - group DR')
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
subplot(2,2,2); imagesc(M_Fsub(dpRSN,dpRSN),[-1 1]); title('Full netmat - subject DR')
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
subplot(2,2,3); imagesc(M_Pgrp(dpRSN,dpRSN),[-1 1]); title('Partial netmat - group DR')
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
subplot(2,2,4); imagesc(M_Psub(dpRSN,dpRSN),[-1 1]); title('Partial netmat - subject DR')
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
set(gcf,'Position',[300 300 900 750],'PaperPositionMode','auto')

% Run CCA on spatial maps from DR on rest data driven by group ICA maps:
subs = dir('Results/DRtask_groupmaps/stage2_*.dtseries.nii'); subs(166) = [];
maps_group_all = zeros(91282,Ncomponents,size(subs,1));
for i = 1:size(subs,1)
    fprintf('Loading map for subject %d\n',i);
    M = ft_read_cifti(fullfile('Results/DRtask_groupmaps',subs(i).name));
    M = M.dtseries; M(isnan(M(:,1))==1,:) = [];
    maps_group_all(:,:,i) = M;
end
C = zeros(size(maps_group_all,3),size(maps_group_all,3));
for n = 1:size(maps_group_all,2)
    fprintf('Calculating spatial input data for dimension %d out of %d\n',n,size(maps_group_all,2));
    conf = demean(conf);
    M = (squeeze(maps_group_all(:,n,:))'-conf*(pinv(conf)*demean(squeeze(maps_group_all(:,n,:)),2)'))';
    C = C + M' * M;
end
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
A = V(:,1:Nkeep); CCAinput.DR_groupmaps = A;
if Nperm>0 ; [~,JBgrotU_group,JBgrotV,JBgrotR,Rmax,Imax,grotRpval,Pmax] = S_CCA(A,Nperm); else [~,JBgrotU_group,JBgrotV,JBgrotR,Rmax,Imax] = S_CCA(A,Nperm); end
%CCAout(5,1) = JBgrotR(1); CCAout(5,2) = Rmax; CCAout(5,3) = grotRpval(1);
maps_group = mean(maps_group_all,3);
maps_group_CCA = zeros(size(maps_group));
for n = 1:Ncomponents; maps_group_CCA(:,n) = corr(squeeze(maps_group_all(:,n,:))',JBgrotU_group(:,1)); end

% Run CCA on spatial maps from DR on rest driven by subject icaW*contrast maps:
maps_subjects_all = zeros(91282,Ncomponents,size(subs,1));
for i = 1:size(subs,1)
    fprintf('Loading map for subject %d\n',i);
    M = ft_read_cifti(fullfile('Results/DRtask_subjectmaps',subs(i).name));
    M = M.dtseries; M(isnan(M(:,1))==1,:) = [];
    maps_subjects_all(:,:,i) = M;
end
C = zeros(size(maps_subjects_all,3),size(maps_subjects_all,3));
for n = 1:size(maps_subjects_all,2)
    fprintf('Calculating spatial input data for dimension %d out of %d\n',n,size(maps_subjects_all,2));
    conf = demean(conf);
    M = (squeeze(maps_subjects_all(:,n,:))'-conf*(pinv(conf)*demean(squeeze(maps_subjects_all(:,n,:)),2)'))';
    C = C + M' * M;
end
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
A = V(:,1:Nkeep); CCAinput.DR_subjectmaps = A;
if Nperm>0 ; [~,JBgrotU_subjects,JBgrotVR,JBgrotR,Rmax,Imax,grotRpval,Pmax] = S_CCA(A,Nperm); else [~,JBgrotU_subjects,JBgrotV,JBgrotR,Rmax,Imax] = S_CCA(A,Nperm); end
%CCAout(7,1) = JBgrotR(1); CCAout(7,2) = Rmax; CCAout(7,3) = grotRpval(1);
maps_subjects = mean(maps_subjects_all,3);
maps_subjects_CCA = zeros(size(maps_subjects));
for n = 1:Ncomponents; maps_subjects_CCA(:,n) = corr(squeeze(maps_subjects_all(:,n,:))',JBgrotU_subjects(:,1)); end

% Run CCA on spatial maps icaW*subject contrast maps:
subs = dir('Results/DRtask_subjectmaps/TaskICAsubject*.dtseries.nii'); subs(166) = [];
maps_Tsubjects_all = zeros(91282,Ncomponents,size(subs,1));
for i = 1:size(subs,1)
    fprintf('Loading map for subject %d\n',i);
    M = ft_read_cifti(fullfile('Results/DRtask_subjectmaps',subs(i).name));
    M = M.dtseries; M(isnan(M(:,1))==1,:) = [];
    maps_Tsubjects_all(:,:,i) = M;
end
C = zeros(size(maps_Tsubjects_all,3),size(maps_Tsubjects_all,3));
for n = 1:size(maps_Tsubjects_all,2)
    fprintf('Calculating spatial input data for dimension %d out of %d\n',n,size(maps_Tsubjects_all,2));
    conf = demean(conf);
    M = (squeeze(maps_Tsubjects_all(:,n,:))'-conf*(pinv(conf)*demean(squeeze(maps_Tsubjects_all(:,n,:)),2)'))';
    C = C + M' * M;
end
[V,D] = eig(C); [~,inds] = sort(diag(D), 'descend'); D = D(inds,inds); V = V(:,inds);
A = V(:,1:Nkeep); CCAinput.icaWtask_subjectmaps = A;
if Nperm>0 ; [~,JBgrotU_Tsubjects,JBgrotVR,JBgrotR,Rmax,Imax,grotRpval,Pmax] = S_CCA(A,Nperm); else [~,JBgrotU_Tsubjects,JBgrotV,JBgrotR,Rmax,Imax] = S_CCA(A,Nperm); end
%CCAout(6,1) = JBgrotR(1); CCAout(6,2) = Rmax; CCAout(6,3) = grotRpval(1);
maps_Tsubjects = mean(maps_Tsubjects_all,3);
maps_Tsubjects_CCA = zeros(size(maps_Tsubjects));
for n = 1:Ncomponents; maps_Tsubjects_CCA(:,n) = corr(squeeze(maps_Tsubjects_all(:,n,:))',JBgrotU_Tsubjects(:,1)); end

% Save CCA inputs
save('Results/CCAinputs_taskStuff.mat','CCAinput')

% Plot map contributions to CCA
figure; plot([sum(abs(maps_group_CCA))' sum(abs(maps_subjects_CCA))' sum(abs(maps_Tsubjects_CCA))']);
legend({'maps group DR rest','maps subjects DR rest','maps icaW*contrasts task'})
set(gca,'xtick',1:15); set(gcf,'Position',[200 200 800 400])
title('Contribution of different maps to CCA results'); ylabel('abs summed correlation across grayordinates')

% Save group-averaged maps for DR stage 2 and for subject maps:
IC = ft_read_cifti('/vols/Data/HCP/Phase2/group900/groupICA/groupICA_3T_HCP820_MSMAll_d15.ica/melodic_IC.dtseries.nii');
IC.dtseries = repmat(IC.dtseries(:,1),1,Ncomponents*3); IC.hdr.dim(6) = Ncomponents*3; IC.time = 1:Ncomponents*3;
I = 1; for n = 1:Ncomponents; IC.dtseries(isnan(IC.dtseries(:,1))==0,I) = maps_group(:,n); I = I+1; IC.dtseries(isnan(IC.dtseries(:,1))==0,I) = maps_subjects(:,n); I = I+1; IC.dtseries(isnan(IC.dtseries(:,1))==0,I) = maps_Tsubjects(:,n); I = I+1; end
ft_write_cifti('Results/DRtaskmaps_averaged',IC,'parameter','dtseries')

% Plot scatter plots:
addpath ~steve/matlab/FACS
figure
subplot(2,2,1); dscatter(Fgrp(:),Fsub(:)); xlabel('Netmats from group maps'); ylabel('Netmats from subject maps'); title('Full netmats - all subject edges'); hold on; plot(-2:0.1:2,-2:0.1:2)
subplot(2,2,2); dscatter(Pgrp(:),Psub(:)); xlabel('Netmats from group maps'); ylabel('Netmats from subject maps'); title('Partial netmats - all subject edges'); hold on; plot(-2:0.1:2,-2:0.1:2)
subplot(2,2,3); scatter(M_Fgrp(:),M_Fsub(:)); xlabel('Netmats from group maps'); ylabel('Netmats from subject maps'); title('Full netmats - mean edges'); hold on; plot(-2:0.1:2,-2:0.1:2)
subplot(2,2,4); scatter(M_Pgrp(:),M_Psub(:)); xlabel('Netmats from group maps'); ylabel('Netmats from subject maps'); title('Partial netmats - mean edges'); hold on; plot(-2:0.1:2,-2:0.1:2)
set(gcf,'Position',[300 300 1000 900],'PaperPositionMode','auto')

% Average maps according to CCA results
D = ft_read_cifti('/vols/Scratch/janineb/HCP/CCA/CCA_new_PFMs/Results/Maps/SynthMap_S820_M50_Aug16_014.dtseries.nii');
D.time = 1:Ncomponents*6; D.hdr.dim(6) = Ncomponents*6; D.dtseries = repmat(D.dtseries(:,1),1,Ncomponents*6);
[~,IU] = sort(JBgrotU_ICA(:,1)); I = 1;
for n = 1:Ncomponents
    G = squeeze(maps_group_all(:,n,IU)); G = G(:,[1:50 501:550]);
    S = squeeze(maps_subjects_all(:,n,IU)); S = S(:,[1:50 501:550]);
    T = squeeze(maps_Tsubjects_all(:,n,IU)); T = T(:,[1:50 501:550]);
    D.dtseries(isnan(D.dtseries(:,1))==0,I) = mean(G(:,1:50),2); I=I+1;
    D.dtseries(isnan(D.dtseries(:,1))==0,I) = mean(G(:,51:100),2); I=I+1;
    D.dtseries(isnan(D.dtseries(:,1))==0,I) = mean(S(:,1:50),2); I=I+1;
    D.dtseries(isnan(D.dtseries(:,1))==0,I) = mean(S(:,51:100),2); I=I+1;
    D.dtseries(isnan(D.dtseries(:,1))==0,I) = mean(T(:,1:50),2); I=I+1;
    D.dtseries(isnan(D.dtseries(:,1))==0,I) = mean(T(:,51:100),2); I=I+1;
end
ft_write_cifti('Results/DRtaskmaps_CCAextremes',D,'parameter','dtseries');

Ms = [mean(JBgrotU_group(IU(501:550),1)) mean(JBgrotU_group(IU(1:50),1)) mean(JBgrotU_subjects(IU(1:50),1)) mean(JBgrotU_subjects(IU(741:790),1))];
Ss = [std(JBgrotU_group(IU(1:50),1)) std(JBgrotU_group(IU(501:550),1)) std(JBgrotU_subjects(IU(1:50),1)) std(JBgrotU_subjects(IU(741:790),1))];
figure; barweb(Ms,Ss);

% Plot CCA results
figure; imagesc(corr([JBgrotU_ICA(:,1:3) JBgrotU_group(:,1:3) JBgrotU_subjects(:,1:3) JBgrotU_Tsubjects(:,1:3)]),[-1 1])
set(gca,'ytick',1:12,'yticklabel',{'THE CCA 1','THE CCA 2','THE CCA3','Group DR 1','Group DR 2','Group DR 3','Subjects DR 1','Subjects DR 2','Subjects DR 3','Subjects icaW 1','Subjects icaW 2','Subjects icaW 3'});
title('Correlations of U''s from different CCA results')
hold on; hline([3.5 6.5 9.5],'k'); vline([3.5 6.5 9.5],'k');

% Correlate netmats and spatial maps against THE CCA result
A = Fgrp; A(166,:) = []; CCA.Fgrp = corr(JBgrotU_ICA(:,1),A);
A = Fsub; A(166,:) = []; CCA.Fsub = corr(JBgrotU_ICA(:,1),A);
A = Pgrp; A(166,:) = []; CCA.Pgrp = corr(JBgrotU_ICA(:,1),A);
A = Psub; A(166,:) = []; CCA.Psub = corr(JBgrotU_ICA(:,1),A);
figure; 
subplot(2,2,1); A = reshape(CCA.Fgrp,Ncomponents,Ncomponents); imagesc(A(dpRSN,dpRSN),[-0.5 0.5]); title('Full netmat group'); ylabel('Correlation with U of THE CCA')
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
subplot(2,2,2); A = reshape(CCA.Fsub,Ncomponents,Ncomponents); imagesc(A(dpRSN,dpRSN),[-0.5 0.5]); title('Full netmat subjects'); ylabel('Correlation with U of THE CCA')
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
subplot(2,2,3); A = reshape(CCA.Pgrp,Ncomponents,Ncomponents); imagesc(A(dpRSN,dpRSN),[-0.5 0.5]); title('Partial netmat group'); ylabel('Correlation with U of THE CCA')
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
subplot(2,2,4); A = reshape(CCA.Psub,Ncomponents,Ncomponents); imagesc(A(dpRSN,dpRSN),[-0.5 0.5]); title('Partial netmat subjects'); ylabel('Correlation with U of THE CCA');
set(gca,'xtick',1:15,'xticklabel',dpRSN,'ytick',1:15,'yticklabel',dpRSN);
set(gcf,'Position',[300 300 900 750],'PaperPositionMode','auto')

CCA.Mgrp = zeros(size(maps_Tsubjects));
for n = 1:Ncomponents; CCA.Mgrp(:,n) = corr(squeeze(maps_group_all(:,n,:))',JBgrotU_ICA(:,1)); end
CCA.Msub = zeros(size(maps_Tsubjects));
for n = 1:Ncomponents; CCA.Msub(:,n) = corr(squeeze(maps_subjects_all(:,n,:))',JBgrotU_ICA(:,1)); end
CCA.Mtsub = zeros(size(maps_Tsubjects));
for n = 1:Ncomponents; CCA.Mtsub(:,n) = corr(squeeze(maps_Tsubjects_all(:,n,:))',JBgrotU_ICA(:,1)); end
IC = ft_read_cifti('/vols/Data/HCP/Phase2/group900/groupICA/groupICA_3T_HCP820_MSMAll_d15.ica/melodic_IC.dtseries.nii');
IC.dtseries = repmat(IC.dtseries(:,1),1,4); IC.hdr.dim(6) = 4; IC.time = 1:4;
IC.dtseries(isnan(IC.dtseries(:,1))==0,1) = mean(abs(CCA.Mgrp),2);
IC.dtseries(isnan(IC.dtseries(:,1))==0,2) = mean(abs(CCA.Msub),2);
IC.dtseries(isnan(IC.dtseries(:,1))==0,3) = mean(abs(CCA.Mtsub),2);
IC.dtseries(isnan(IC.dtseries(:,1))==0,4) = mean(abs([CCA.Mgrp CCA.Msub CCA.Mtsub]),2);
ft_write_cifti('Results/DRtaskmaps_CCAaveraged',IC,'parameter','dtseries')

% Summary of CCA correlations
ms = [mean(abs(CCA.Mgrp(:))) mean(abs(CCA.Msub(:))) mean(abs(CCA.Mtsub(:))); ...
    nanmean(abs(CCA.Fgrp(:))) nanmean(abs(CCA.Fsub(:))) nan; ...
    nanmean(abs(CCA.Pgrp(:))) nanmean(abs(CCA.Psub(:))) nan];
ss = [std(abs(CCA.Mgrp(:))) std(abs(CCA.Msub(:))) std(abs(CCA.Mtsub(:))); ...
    nanstd(abs(CCA.Fgrp(:))) nanstd(abs(CCA.Fsub(:))) nan; ...
    nanstd(abs(CCA.Pgrp(:))) nanstd(abs(CCA.Psub(:))) nan];
figure; barweb(ms,ss)





