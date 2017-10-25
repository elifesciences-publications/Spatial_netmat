% This script gets called by CCArun_all_inputs_RUN.sh

% Set and load inputs:
Dir = '/home/fs0/janineb/scratch/';
Ddir = Ddir(1:end-4);
S = fullfile(Dir,'HCP','netmat_simulations');
%Nperm = 100000;
Nperm = 0;

% Set up paths
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;

brainsAll = load(fullfile(S,'Results',Ddir));
A1 = fieldnames(brainsAll);
OUTPUT_CCA = zeros(size(A1,1),6);
for x = 1:size(A1,1)
    fprintf('running CCA on %s.%s\n',Ddir,A1{x});
    brains = brainsAll.(A1{x});
    
    % Run CCA:
    if Nperm>0
        [~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.Rmax,CCA.Imax,CCA.grotRpval,CCA.Pmax] = S_CCA(brains,Nperm);
    else
        [~,CCA.JBgrotU,CCA.JBgrotV,CCA.JBgrotR,CCA.Rmax,CCA.Imax] = S_CCA(brains,Nperm);
        CCA.grotRpval = zeros(size(CCA.JBgrotR)); CCA.Pmax = 0;
    end
    
    OUTPUT_CCA(x,:) = [CCA.JBgrotR(1) CCA.Rmax CCA.Imax CCA.Pmax CCA.grotRpval(1) sum(CCA.grotRpval<0.05)];
    
end

% Load existing output and add CCA results to it
load(fullfile(S,'Results',sprintf('results_%s.mat',Ddir(11:end))),'OUTPUT','F','P');
save(fullfile(S,'Results',sprintf('results_%s.mat',Ddir(11:end))),'OUTPUT','OUTPUT_CCA','F','P','-v7.3');




