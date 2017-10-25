%clear all; close all; clc

FMRIB='/vols/Data/HCP/Phase2/scripts900';
subs_all = dir('/home/fs0/samh/scratch/HCP/S820_M50_Aug16.pfm/FinalModel/Subjects');
subs_all = subs_all(3:end); subs_all(169) = [];
S = zeros(size(subs_all,1),2); for n = 1:size(subs_all,1); S(n,1) = str2double(subs_all(n).name); end; subs_all = S; clear n
%load('Results/subs_Matt.mat')
%EB = hcp2blocks(fullfile(FMRIB,'vars','RESTRICTED_fmribsteve_12_14_2015_10_10_6.csv'),'',false,subs_all);

ignore_sign = 0;
AllBoots = zeros(size(EB,1),1000);
for X = 1:1000
    
    % First level (pick family type):
    if ignore_sign == 1; repl = true; else repl = false; end
    level1 = randsample(EB(:,2),length(EB),repl);
    
    % Second level (pick specific family):
    N = unique(level1);
    level2 = [];
    for n = N'
        if ignore_sign==1 || sign(n)==1;
            repl=true;
        elseif sign(n)==-1;
            repl=false;
        end
        A = EB(EB(:,2)==n,3); if length(A)==1; A = [A; A]; end
        level2 = [level2; randsample(A,length(find(level1==n)),true)];
    end
    
    % Third level (pick sibling/ monozygote/ dizygote pair):
    N = unique(level2);
    level3 = [];
    for n = N'
        %%% NOT POSSIBLE TO SAMPLE WITHOUT REPLACEMENT BECAUSE OFTEN MORE
        %%% SUBJECTS THAN IN FAMILY (DUE TO REPLACEMENT AT PREVIOUS LEVEL)
        %     if ignore_sign==1 || sign(n)==1;
        %         repl=true;
        %     elseif sign(n)==-1;
        %         repl=false;
        %     end
        repl = true;
        A = EB(EB(:,3)==n,4); if length(A)==1; A = [A; A]; end
        level3 = [level3; randsample(A,length(find(level2==n)),repl)];
    end
    
    % Fourth level (pick specific subjects):
    N = unique(level3);
    level4 = [];
    for n = N'
        if ignore_sign==1 || sign(n)==1;
            repl=true;
        elseif sign(n)==-1;
            repl=false;
        end
        A = EB(EB(:,4)==n,5); if length(A)==1; A = [A; A]; end
        level4 = [level4; randsample(A,length(find(level3==n)),repl)];
    end
    for n = 1:size(S,1)
        AllBoots(n,X) = find(S(:,1)==level4(n));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set and load inputs:
Ddir = 'input_CCA_table1.mat';
Dir = '/home/fs0/janineb/scratch/';
Ddir = Ddir(1:end-4);
Nperm = 0;

% Set up paths
path(path,fullfile(Dir,'matlab'))
path(path,fullfile(Dir,'HCP','CCA','files'))
path(path,fullfile(Dir,'matlab','cifti-matlab'))
path(path,fullfile(Dir,'HCP','DMN','DMN_functions'))
addpath ~steve/NETWORKS/FSLNets;
path(path,fullfile(Dir,'HCP','netmat_simulations','palm-alpha106'));
S = fullfile('/vols/Scratch/janineb/','HCP','CCA');

brainsAll = load(fullfile('Results',Ddir));
A1 = fieldnames(brainsAll);
R = zeros(1000,size(A1,1));
for x = 1:size(A1,1)
    fprintf('running CCA on %s.%s\n',Ddir,A1{x});
    brains = brainsAll.(A1{x});
    if size(brains,1)==819
        load(fullfile(S,'files','Permutation_100000.mat'),'uu2');
        for n = 1:1000
            [~,~,r,~,~,~] = canoncorr(brains(AllBoots(:,n),:),uu2(AllBoots(:,n),:));
            R(n,x) = r(1);
            [~,~,r,~,~,~] = canoncorr(brains,uu2);
        end
    end
end

