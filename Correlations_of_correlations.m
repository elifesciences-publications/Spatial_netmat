clear all; close all; clc

%%% Do we need to estimate the inputs or just plot?
estimate = 0;
plot = 1;

% Set up paths
Dir = '/home/fs0/janineb/scratch/';
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;
Sdir = '/vols/Scratch/janineb/HCP/CCA/';
HCPdir = '/vols/Data/HCP/Phase2/group900';
subs = dir('/vols/Data/HCP/Phase2/subjects900/*'); subs = subs(3:end);

if estimate == 1;
    
    %%% Get averaged U to find shared variance
    brainsAll = load('Results/input_CCA_table1.mat');
    A1 = fieldnames(brainsAll);
    U = []; Rmax = []; names = [];
    for x = 1:size(A1,1)
        brains = brainsAll.(A1{x});
        [~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.Rmax,CCA.Imax] = S_CCA(brains,0);
        U = [U CCA.JBgrotU(:,CCA.Imax)];
        Rmax = [Rmax; CCA.Rmax CCA.Imax];
        names = [names; A1(x)];
    end
    O = {'025','200','Yeo'};
    for x = 1:length(O)
        brainsAll = load(sprintf('Results/input_CCA_%s.mat',O{x}));
        A1 = fieldnames(brainsAll);
        for n = [2 4]
            brains = brainsAll.(A1{n});
            [~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.Rmax,CCA.Imax] = S_CCA(brains,0);
            U = [U CCA.JBgrotU(:,CCA.Imax)];
            Rmax = [Rmax; CCA.Rmax CCA.Imax];
            names = [names; [O{x} A1{n}]];
        end
    end
    ingrot = U;
    ingrot = demean(ingrot); CORR.U = corr(ingrot');
    
    %%% Get subject correlation matrix across all behaviour
    U = load('~samh/DataSets/HCP/BehaviouralAnalyses/ImputedSoftImpute_Variables.txt');
    N = readtable('/home/fs0/samh/DataSets/HCP/BehaviouralAnalyses/ImputedSoftImpute_VariableHeaders.txt');
    N = table2cell(N); I = strncmp('FS',N,2); I = find(I==1)+1;
    U(169,:) = [];
    U1 = demean(U(:,I)); U2 = demean(U(:,setdiff(1:400,I)));
    CORR.behaviour = corr(U2'); CORR.FS = corr(U1');
    
    %%% Input for timeseries, netmat and amplitudes
    in = 'MSMAll';
    d = 200;
    Ddir = sprintf('3T_HCP820_%s_d%d_ts2',in,d);
    
    %%% Timeseries
    fprintf('Loading ICA200 timeseries \n');
    ts = nets_load(fullfile(HCPdir,'node_timeseries',Ddir),0.72,0,1);
    
    %%% Netmats
    fprintf('Calculating ICA200 netmats \n');
    Pnetmat = nets_netmats(ts,1,'ridgep',0.01);
    ingrot = Pnetmat;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA200_Pnetmat = corr(ingrot');
    Fnetmat = nets_netmats(ts,1,'corr');
    ingrot = Fnetmat;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA200_Fnetmat = corr(ingrot');
    clear Pnetmat Fnetmat ingrot
    
    %%% Amplitudes
    fprintf('Calculating ICA200 amplitudes \n');
    [~,all_stats] = nets_stats(ts);
    ingrot = all_stats.std;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA200_amp = corr(ingrot');
    clear ingrot
    
    %%% Shifts ICA DR maps
    fprintf('Calculating ICA200 spatial maps \n');
    maps = zeros(size(subs,1),91282,d);
    for i = 1:size(subs,1)
        M = ft_read_cifti(fullfile(HCPdir,'node_maps',Ddir,[subs(i).name '.dtseries.nii']));
        M = M.dtseries; M(isnan(M(:,1))==1,:) = [];
        maps(i,:,:) = M;
    end
    for i = 1:d
        m = maps(:,:,i)';
        m = demean(m,2);
        s = demean(all_stats.std(:,i));
        m = (m'-s*(pinv(s)*m'))';
        maps(:,:,i) = m';
    end
    maps = reshape(maps,size(maps,1),size(maps,2)*size(maps,3));
    ingrot = maps;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA200_spatial = corr(ingrot');
    clear maps ingrot M s all_stats i ts
    
    %%% ICA 25 INPUTS
    d = 25;
    Ddir = sprintf('3T_HCP820_%s_d%d_ts2',in,d);
    ts = nets_load(fullfile(HCPdir,'node_timeseries',Ddir),0.72,0,1);
    
    %%% Netmats d=25
    fprintf('Calculating ICA25 netmats \n');
    Pnetmat = nets_netmats(ts,1,'ridgep',0.01);
    ingrot = Pnetmat;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA25_Pnetmat = corr(ingrot');
    Fnetmat = nets_netmats(ts,1,'corr');
    ingrot = Fnetmat;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA25_Fnetmat = corr(ingrot');
    clear Pnetmat Fnetmat ingrot
    
    %%% Amplitudes d=25
    fprintf('Calculating ICA25 amplitudes \n');
    [~,all_stats] = nets_stats(ts);
    ingrot = all_stats.std;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA25_amp = corr(ingrot');
    clear ingrot
    
    %%% Shifts ICA DR maps d=25
    fprintf('Calculating ICA5 spatial maps \n');
    maps = zeros(size(subs,1),91282,d);
    for i = 1:size(subs,1)
        M = ft_read_cifti(fullfile(HCPdir,'node_maps',Ddir,[subs(i).name '.dtseries.nii']));
        M = M.dtseries; M(isnan(M(:,1))==1,:) = [];
        maps(i,:,:) = M;
    end
    for i = 1:d
        m = maps(:,:,i)';
        m = demean(m,2);
        s = demean(all_stats.std(:,i));
        m = (m'-s*(pinv(s)*m'))';
        maps(:,:,i) = m';
    end
    maps = reshape(maps,size(maps,1),size(maps,2)*size(maps,3));
    ingrot = maps;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.ICA25_spatial = corr(ingrot');
    clear maps ingrot M s all_stats i ts
    
    %%% Shifts (PFM subject maps)
    fprintf('Loading PFM data \n');
    Ddir = 'S820_M50_Aug16.pfm';
    d = 50;
    PFMdir = '/home/fs0/samh/scratch/HCP/';
    P = fullfile(PFMdir,Ddir,'FinalModel');
    noise = [5 15 16 21 25 28 32 36 41 44];
    signal = setdiff(1:d,noise);
    maps = PFM_loadSubjectSpatialMaps(P,signal); % 820 x 91282 x d
    maps = reshape(maps,size(maps,1),size(maps,2)*size(maps,3));
    ingrot = maps;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.PFM50_spatial = corr(ingrot');
    clear maps ingrot
    
    %%% PFM amplitudes
    fprintf('Calculating PFM50 amplitudes \n');
    runs = {'LR','RL'};
    subs = dir(fullfile(P,'Subjects'));
    A = zeros(820,d*4);
    for i = 1:size(subs,1)-2
        I = 1;
        for day = 1:2
            for r = 1:2
                A(i,(I-1)*d+1:I*d) = h5read(fullfile(P,'Subjects',subs(i+2).name,'Runs',sprintf('%d_%s',day,runs{r}),'ComponentWeightings.post','Means.hdf5'),'/dataset');
                I = I+1;
            end
        end
    end
    ingrot = A;
    ingrot(169,:) = []; ingrot = demean(ingrot); CORR.PFM50_amp = corr(ingrot');
    clear runs subs A i I day r ingrot
    
    %%% Save output
    save('Results/Correlations_of_correlations.mat','CORR','-v7.3');
end

if plot == 1
    load('Results/Correlations_of_correlations.mat')
    AllIn = []; Names = [];
    A1 = fieldnames(CORR);
    Inets = [2 3 6 7];
    Iother = [9 5 10 11];
    figure; set(gcf,'position',[20 20 908 767],'PaperPositionMode','auto')
    for x = 1:length(Inets)
        a = CORR.(A1{Inets(x)});
        a = a(tril(ones(size(a)),-1)==1);
        AllIn =  0.5*log((1+a)./(1-a));
        Names = A1(Inets(x));
        for n = 1:length(Iother)
            a = CORR.(A1{Iother(n)});
            a = a(tril(ones(size(a)),-1)==1);
            AllIn = [AllIn  0.5*log((1+a)./(1-a));];
            Names = [Names; A1(Iother(n))];
        end
        Names = strrep(Names,'_',' '); Names = strrep(Names,'amp','amplitude'); 
        Names = strrep(Names,'Pnetmat','partial netmat'); Names = strrep(Names,'Fnetmat','full netmat');
        Names = strrep(Names,'netmat','network matrix');
        N = size(AllIn,2);
        P = cov(AllIn);
        rho = 0.01; P = P/sqrt(mean(diag(P).^2)); P = inv(P+rho*eye(N)); P = -P;
        P = (P./repmat(sqrt(abs(diag(P))),1,N)) ./ repmat(sqrt(abs(diag(P)))',N,1);
        F = corr(AllIn);
        OUT = tril(F)+triu(P);
        OUT(eye(N)>0) = inf; colormap('jet'); grotc=colormap; grotc(end,:)=[.8 .8 .8]; colormap(grotc);
        subplot(2,2,x); imagesc(OUT,[-1.1 1.1]);
        set(gca,'ytick',1:N,'yticklabel',Names,'xtick',1:N,'xticklabel',Names)
        %XYrotalabel(20);
        title(sprintf('%s',Names{1}));
        rectangle('position',[1.5 0.52 4 0.98],'EdgeColor','k','LineWidth',2);
        rectangle('position',[0.52 1.5 0.98 4],'EdgeColor','k','LineWidth',2);
    end
    print(gcf,'-dtiff','-r300','Results/Correlation_correlation_matrices')
    colorbar('location','southoutside');
    print(gcf,'-dtiff','-r300','Results/Correlation_correlation_matrices_colorbar')
    
    load('Results/results_200.mat','F');
    A = reshape(F.netORIG(1,:),200,200); A = tril(A,-1); A(A==0) = inf; I = find(isinf(A)==0);
    figure; imagesc(A,[-1 1]); xlabel('Nodes'); ylabel('Nodes'); colormap(grotc); colorbar; set(gca,'xtick',50:50:200,'ytick',50:50:200);
    print(gcf,'-dtiff','-r150','Results/Correlation_correlation_matrices_explanation1')
    figure; imagesc(F.netORIG(:,I),[-1 1]); xlabel('Edges'); ylabel('Subjects'); colorbar; set(gca,'xtick',2000:4000:19900,'ytick',200:200:800);
    print(gcf,'-dtiff','-r150','Results/Correlation_correlation_matrices_explanation2')
    A = F.netORIG(:,I); B = F.netNEW(:,I); A = corr(A'); A = tril(A,-1); A(A==0) = inf; I = find(isinf(A)==0);
    figure; imagesc(A,[0 1]); xlabel('Subjects'); ylabel('Subjects'); colormap(grotc); colorbar; set(gca,'xtick',200:200:800,'ytick',200:200:800);
    print(gcf,'-dtiff','-r150','Results/Correlation_correlation_matrices_explanation3')
    A = A(I); B = corr(B'); B = B(I);
    figure; imagesc(A,[0 1]); axis off; set(gcf,'Position',[0 0 50 1000],'PaperPositionMode','auto')
    print(gcf,'-dtiff','-r150','Results/Correlation_correlation_matrices_explanation4')
    figure; imagesc(B,[0 1]); axis off; set(gcf,'Position',[0 0 50 1000],'PaperPositionMode','auto')
    print(gcf,'-dtiff','-r150','Results/Correlation_correlation_matrices_explanation5')
end




