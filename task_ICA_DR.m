
SUBJECTS = '/vols/Scratch/HCP/rfMRI/subjects900';
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
addpath ~steve/NETWORKS/FSLNets;
P = '/vols/Scratch/HCP/tfMRI/Q900/';
MSM = '_MSMAll';
Ncomponents = 15;

ff{1}=sprintf('%s/%s/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii',SUBJECTS,subID);
ff{2}=sprintf('%s/%s/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii',SUBJECTS,subID);
ff{3}=sprintf('%s/%s/MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii',SUBJECTS,subID);
ff{4}=sprintf('%s/%s/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii',SUBJECTS,subID);

tt = {'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};

% Load ICA results
load(sprintf('Results/task_ICA_%d.mat',Ncomponents),'icaW');

% Load subject data
for i=1:4
  BO=ft_read_cifti(ff{i});  
  BO = BO.dtseries; BO(isnan(BO(:,1))==1,:) = [];
  SUB{i}=BO; clear BO;
end

% Load subject map and perform DR
if exist(fullfile(P,subID,'MNINonLinear','Results',...
        sprintf('%s_All_task_maps_%s.dtseries.nii',MSM(2:end),subID)),'file');
    BO = ft_read_cifti(fullfile(P,subID,'MNINonLinear','Results',...
        sprintf('%s_All_task_maps_%s.dtseries.nii',MSM(2:end),subID)));
    BO = BO.dtseries; BO(isnan(BO(:,1))==1,:) = [];
    BO = (icaW * BO')';
    IC = ft_read_cifti('/vols/Data/HCP/Phase2/group900/groupICA/groupICA_3T_HCP820_MSMAll_d15.ica/melodic_IC.dtseries.nii');
    IC.dtseries = repmat(IC.dtseries(:,1),1,size(BO,2)); IC.hdr.dim(6) = size(BO,2); IC.time = 1:size(BO,2);
    IC.dtseries(isnan(IC.dtseries(:,1))==0,:) = BO;
    ft_write_cifti(sprintf('Results/DRtask_subjectmaps/TaskICAsubject_%s',subID),IC,'parameter','dtseries')
    BO = BO./repmat(max(BO),size(BO,1),1); % Peak norm maps before feeding into DR
    pGM = pinv(nets_demean(BO));  allNODEts = []; allNODEts2 = []; allMAPS=zeros(size(BO));
    clear D N n IC;
    for i = 1:4
        NODEts = nets_demean((pGM*nets_demean(SUB{i}))');
        allNODEts = [allNODEts ; NODEts];
        allNODEts2 = [allNODEts2; nets_demean((BO'*SUB{i})')];
        NODEts = nets_normalise(NODEts); pTSd{i}=pinv(NODEts)'; grot=nets_demean(SUB{i}');
        allMAPS = allMAPS + ( grot' * pTSd{i} ) ;
    end
    clear grot pTSd;
    dlmwrite(sprintf('Results/DRtask_subjectmaps/%s.txt',subID),allNODEts,'delimiter',' ');
    dlmwrite(sprintf('Results/DRtask_subjectmaps/masked_ts/%s_masked.txt',subID),allNODEts2,'delimiter',' ');
    clear allNODEts NODEts allNODEts2;
    IC = ft_read_cifti('/vols/Data/HCP/Phase2/group900/groupICA/groupICA_3T_HCP820_MSMAll_d15.ica/melodic_IC.dtseries.nii');
    IC.dtseries = repmat(IC.dtseries(:,1),1,size(allMAPS,2)); IC.hdr.dim(6) = size(allMAPS,2); IC.time = 1:size(allMAPS,2);
    IC.dtseries(isnan(IC.dtseries(:,1))==0,:) = allMAPS/4;
    ft_write_cifti(sprintf('Results/DRtask_subjectmaps/stage2_%s',subID),IC,'parameter','dtseries')
    clear allMAPS;   
    
    % Load group map and perform DR
    BO = ft_read_cifti('/vols/Scratch/janineb/HCP/netmat_simulations/Task_ICA15.dtseries.nii');
    BO = BO.dtseries; BO(isnan(BO(:,1))==1,:) = [];
    BO = BO./repmat(max(BO),size(BO,1),1); % Peak norm maps before feeding into DR
    pGM = pinv(nets_demean(BO));  allNODEts = []; allMAPS=zeros(size(BO));
    clear D N n IC;
    for i = 1:4
        NODEts = nets_demean((pGM*nets_demean(SUB{i}))');
        allNODEts = [allNODEts ; NODEts];
        NODEts=nets_normalise(NODEts); pTSd{i}=pinv(NODEts)'; grot=nets_demean(SUB{i}');
        allMAPS = allMAPS + ( grot' * pTSd{i} ) ;
    end
    clear grot pTSd;
    dlmwrite(sprintf('Results/DRtask_groupmaps/%s.txt',subID),allNODEts,'delimiter',' ');
    clear allNODEts NODEts;
    IC = ft_read_cifti('/vols/Data/HCP/Phase2/group900/groupICA/groupICA_3T_HCP820_MSMAll_d15.ica/melodic_IC.dtseries.nii');
    IC.dtseries = repmat(IC.dtseries(:,1),1,size(allMAPS,2)); IC.hdr.dim(6) = size(allMAPS,2); IC.time = 1:size(allMAPS,2);
    IC.dtseries(isnan(IC.dtseries(:,1))==0,:) = allMAPS/4;
    ft_write_cifti(sprintf('Results/DRtask_groupmaps/stage2_%s',subID),IC,'parameter','dtseries')
    clear allMAPS;
    
    % DR on task data:
    allMAPS = zeros(size(BO));
    rr = {'RL','LR'};
    if exist(fullfile(P,'..',subID),'dir');
        for n = 1:length(tt)
            for r = 1:2
                %D = ft_read_cifti(fullfile('/vols/Scratch/janineb/',subID,'MNINonLinear','Results',sprintf('tfMRI_%s_%s',tt{n},rr{r}),sprintf('tfMRI_%s_%s_Atlas_MSMAll.dtseries.nii',tt{n},rr{r})));
                D = ft_read_cifti(fullfile(P,'..',subID,'MNINonLinear','Results',sprintf('tfMRI_%s_%s',tt{n},rr{r}),sprintf('tfMRI_%s_%s_Atlas.dtseries.nii',tt{n},rr{r})));
                D = D.dtseries; D(isnan(D(:,1))==1,:) = [];             
                NODEts = nets_demean((pGM*nets_demean(D))');
                NODEts = nets_normalise(NODEts); pT=pinv(NODEts)'; grot=nets_demean(D');
                allMAPS = allMAPS + ( grot' * pT ) ;
            end
        end
    end
    IC = ft_read_cifti('/vols/Data/HCP/Phase2/group900/groupICA/groupICA_3T_HCP820_MSMAll_d15.ica/melodic_IC.dtseries.nii');
    IC.dtseries = repmat(IC.dtseries(:,1),1,size(allMAPS,2)); IC.hdr.dim(6) = size(allMAPS,2); IC.time = 1:size(allMAPS,2);
    IC.dtseries(isnan(IC.dtseries(:,1))==0,:) = allMAPS/14;
    ft_write_cifti(sprintf('Results/DRtask_groupmaps/stage2task_%s',subID),IC,'parameter','dtseries')
    %clear allMAPS NODEts D pGM; 
end




