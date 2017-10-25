clear all; close all; clc

% Set and load inputs:
Ddir = 'input_CCA_table1_withWarp.mat';
Dir = '/home/fs0/janineb/scratch/';
Ddir = Ddir(1:end-4);
S = fullfile(Dir,'HCP','netmat_simulations');
Nperm = 100000;
%Nperm = 0;

% Load Steve's reference result
U_ICA = load(fullfile('/vols/Scratch/janineb/','HCP','CCA','files','ICA200_MSMall_PartialNetmat_JBgrotU.txt'));

% Set up paths
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;
path(path,fullfile(Dir,'HCP','netmat_simulations','palm-alpha106'));

brainsAll = load(fullfile(S,'Results',Ddir));
A1 = fieldnames(brainsAll);
OUTPUT_CCA = zeros(size(A1,1),6);
thresh05 = zeros(size(A1,1),4);
nulls = struct; Us = struct;
for x = 14%:size(A1,1)
    fprintf('running CCA on %s.%s\n',Ddir,A1{x});
    brains = brainsAll.(A1{x});
    
    % Run CCA:
    if Nperm>0
        %[~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.JBgrotF,CCA.Rmax,CCA.Imax,CCA.grotRpval,CCA.grotFpval,CCA.FPmax,CCA.Fnull] = S_CCA_Fstat(brains,Nperm);
        [~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.Rmax,CCA.Imax,CCA.grotRpval,CCA.Pmax,CCA.null] = S_CCA(brains,Nperm);
%         Gdist = 0.5*log((1+CCA.null(:,1))./(1-CCA.null(:,1)));
%         G = 0.5*log((1+CCA.JBgrotR(1))./(1-CCA.JBgrotR(1)));
%         p = palm_datapval(Gdist,Gdist,false);
%         pareto = palm_pareto(Gdist,Gdist,false,0.1,false);
%         lp = -log10(p); lpareto = -log10(pareto);
%         figure; 
%         subplot(2,2,1); qqplot(lp,lpareto); xlabel('-log10(p)'); ylabel('-log10(pareto)'); title('QQ plot')
%         subplot(2,2,2); plot(sort(Gdist),sort(lp)); hold on; plot(sort(Gdist),sort(lpareto),'r'); xlabel('CCA Fstat'); ylabel('-log10(pareto'); title('tail fit based on log p');
%         subplot(2,2,3); plot(sort(Gdist),sort(1-p)); hold on; plot(sort(Gdist),sort(1-pareto),'r'); xlabel('CCA Fstat'); ylabel('1-pareto'); title('tail fit based on 1-p');
%         subplot(2,2,4); plot(sort(Gdist),sort(1-p)); hold on; plot(sort(Gdist),sort(1-pareto),'r'); xlabel('CCA Fstat'); ylabel('1-pareto'); title('tail fit based on 1-p (zoomed in)'); axis([min(Gdist) max(Gdist) 0.96 1.001])
%         print(gcf,'-dtiff','-r300',sprintf('Results/Pareto_tail_%s.tif',A1{x}))
    else
        [~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.Rmax,CCA.Imax] = S_CCA(brains,Nperm);
        CCA.grotRpval = zeros(size(CCA.JBgrotR)); CCA.Pmax = 0;
    end
    
    % Run CCA against Steve's result to get variance explained:
    A = U_ICA;
    if size(brains,1)==441
        A = [A(1:168,:); ones(1,size(A,2)); A(169:end,:)];
        load('Results/Data_Matt/ts_real_sim_subs.mat','subs_all');
        subs_remove = setdiff(1:820,subs_all);
        A(unique([subs_remove 169]),:) = [];
    elseif size(brains,1)<790
        A = [A(1:168,:); ones(1,size(A,2)); A(169:end,:)];
        subs_remove = [122 160 162 169 176 200 223 248 250 260 286 292 295 312 315 320 337 353 367 425 457 463 607 669 679 684 690 691 701 760 764 797];
        A(subs_remove,:) = [];
    elseif size(brains,1)<800
        A = [A(1:168,:); ones(1,size(A,2)); A(169:end,:)];
        subs_remove = [122 160 162 169 200 248 250 260 286 292 295 312 315 320 337 353 367 425 457 463 607 669 679 684 690 691 701 760 764 797];
        A(subs_remove,:) = [];
    end
    [~,~,R,~,~,~] = canoncorr(CCA.JBgrotU(:,1:3),A(:,1));
    OUTPUT_CCA(x,:) = [CCA.JBgrotR(1) CCA.grotRpval(1) CCA.Rmax sum(CCA.grotRpval<0.05) CCA.Imax CCA.Pmax];
    thresh05(x,:) = prctile(CCA.null(:,1),[90 95 99 99.8]);
    nulls.(A1{x}) = CCA.null(:,1);
    Us.(A1{x}) = CCA.JBgrotU;
end

% Load existing output and add CCA results to it
save(fullfile(S,'Results','results_CCA_table1_NEWER.mat'),'OUTPUT_CCA','A1');
save(fullfile(S,'Results','results_CCA_table1_thresh05.mat'),'thresh05','A1','nulls','Us','-v7.3');



