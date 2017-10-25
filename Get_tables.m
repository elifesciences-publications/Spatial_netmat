clear all; close all; clc

O = {'200','025','Yeo','Glasser','Glasser_group_parcellation_'};
S = '/vols/Scratch/janineb/HCP/netmat_simulations/Results';

Headers = {'Rnetmat','Rcorrelation','CCA_Ruv','CCA_Puv1','CCA_Ruuica_max','CCA_Ruuica_N','CCA_Ruuica_p'};

%%%%% GET TABLE 2:

T2_partial = nan(25,7);
I = 1;
for n = 1:length(O);
    
    % First line for group netmats, amps and spatial maps (GGG)
    A = fullfile(S,sprintf('results_%s_using_group_maps_using_group_amps.mat',O{n}));
    if exist(A,'file');
        load(A)
        T2_partial(I,1) = OUTPUT(2,2);
        T2_partial(I,2) = OUTPUT(2,1);
        if exist('OUTPUT_CCA','var')
            if size(OUTPUT_CCA,1)>2
            T2_partial(I,3) = OUTPUT_CCA(3,1);
            T2_partial(I,4) = OUTPUT_CCA(3,5);
            T2_partial(I,5) = OUTPUT_CCA(3,2);
            T2_partial(I,6) = OUTPUT_CCA(3,3);
            T2_partial(I,7) = OUTPUT_CCA(3,4);
            end
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Second line for group netmats and subject amps & spatial maps (G--)
    A = fullfile(S,sprintf('results_%s.mat',O{n}));
    if exist(A,'file');
        load(A)
        T2_partial(I,1) = OUTPUT(2,2);
        T2_partial(I,2) = OUTPUT(2,1);
        if exist('OUTPUT_CCA','var')
            if size(OUTPUT_CCA,1)>2
            T2_partial(I,3) = OUTPUT_CCA(3,1);
            T2_partial(I,4) = OUTPUT_CCA(3,5);
            T2_partial(I,5) = OUTPUT_CCA(3,2);
            T2_partial(I,6) = OUTPUT_CCA(3,3);
            T2_partial(I,7) = OUTPUT_CCA(3,4);
            end
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Third line for subject netmats and group amps & spatial maps (-GG)
    A = fullfile(S,sprintf('results_%s_using_group_maps_using_group_amps_using_subject_netmats.mat',O{n}));
    if exist(A,'file');
        load(A)
        T2_partial(I,1) = OUTPUT(2,2);
        T2_partial(I,2) = OUTPUT(2,1);
        if exist('OUTPUT_CCA','var')
            if size(OUTPUT_CCA,1)>2
            T2_partial(I,3) = OUTPUT_CCA(3,1);
            T2_partial(I,4) = OUTPUT_CCA(3,5);
            T2_partial(I,5) = OUTPUT_CCA(3,2);
            T2_partial(I,6) = OUTPUT_CCA(3,3);
            T2_partial(I,7) = OUTPUT_CCA(3,4);
            end
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Fourth line for group netmats, subject amps, group spatial maps (G-G)
    A = fullfile(S,sprintf('results_%s_using_group_maps.mat',O{n}));
    if exist(A,'file');
        load(A)
        T2_partial(I,1) = OUTPUT(2,2);
        T2_partial(I,2) = OUTPUT(2,1);
        if exist('OUTPUT_CCA','var')
            if size(OUTPUT_CCA,1)>2
            T2_partial(I,3) = OUTPUT_CCA(3,1);
            T2_partial(I,4) = OUTPUT_CCA(3,5);
            T2_partial(I,5) = OUTPUT_CCA(3,2);
            T2_partial(I,6) = OUTPUT_CCA(3,3);
            T2_partial(I,7) = OUTPUT_CCA(3,4);
            end
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Fifth line for group netmats and amps, and subject spatial maps (GG-)
    A = fullfile(S,sprintf('results_%s_using_group_amps.mat',O{n}));
    if exist(A,'file');
        load(A)
        T2_partial(I,1) = OUTPUT(2,2);
        T2_partial(I,2) = OUTPUT(2,1);
        if exist('OUTPUT_CCA','var')
            if size(OUTPUT_CCA,1)>2
                T2_partial(I,3) = OUTPUT_CCA(3,1);
                T2_partial(I,4) = OUTPUT_CCA(3,5);
                T2_partial(I,5) = OUTPUT_CCA(3,2);
                T2_partial(I,6) = OUTPUT_CCA(3,3);
                T2_partial(I,7) = OUTPUT_CCA(3,4);
            end
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
end

%%%%% GET SUPPLEMENTARY TABLE 1

TS1_full = nan(25,7);
I = 1;
for n = 1:length(O);
    
    % First line for group netmats, amps and spatial maps (GGG)
    A = fullfile(S,sprintf('results_%s_using_group_maps_using_group_amps.mat',O{n}));
    if exist(A,'file');
        load(A)
        TS1_full(I,1) = OUTPUT(1,2);
        TS1_full(I,2) = OUTPUT(1,1);
        if exist('OUTPUT_CCA','var')
            TS1_full(I,3) = OUTPUT_CCA(1,1);
            TS1_full(I,4) = OUTPUT_CCA(1,5);
            TS1_full(I,5) = OUTPUT_CCA(1,2);
            TS1_full(I,6) = OUTPUT_CCA(1,3);
            TS1_full(I,7) = OUTPUT_CCA(1,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Second line for group netmats and subject amps & spatial maps (G--)
    A = fullfile(S,sprintf('results_%s.mat',O{n}));
    if exist(A,'file');
        load(A)
        TS1_full(I,1) = OUTPUT(1,2);
        TS1_full(I,2) = OUTPUT(1,1);
        if exist('OUTPUT_CCA','var')
            TS1_full(I,3) = OUTPUT_CCA(1,1);
            TS1_full(I,4) = OUTPUT_CCA(1,5);
            TS1_full(I,5) = OUTPUT_CCA(1,2);
            TS1_full(I,6) = OUTPUT_CCA(1,3);
            TS1_full(I,7) = OUTPUT_CCA(1,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Third line for subject netmats and group amps & spatial maps (-GG)
    A = fullfile(S,sprintf('results_%s_using_group_maps_using_group_amps_using_subject_netmats.mat',O{n}));
    if exist(A,'file');
        load(A)
        TS1_full(I,1) = OUTPUT(1,2);
        TS1_full(I,2) = OUTPUT(1,1);
        if exist('OUTPUT_CCA','var')
            TS1_full(I,3) = OUTPUT_CCA(1,1);
            TS1_full(I,4) = OUTPUT_CCA(1,5);
            TS1_full(I,5) = OUTPUT_CCA(1,2);
            TS1_full(I,6) = OUTPUT_CCA(1,3);
            TS1_full(I,7) = OUTPUT_CCA(1,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Fourth line for group netmats, subject amps, group spatial maps (G-G)
    A = fullfile(S,sprintf('results_%s_using_group_maps.mat',O{n}));
    if exist(A,'file');
        load(A)
        TS1_full(I,1) = OUTPUT(1,2);
        TS1_full(I,2) = OUTPUT(1,1);
        if exist('OUTPUT_CCA','var')
            TS1_full(I,3) = OUTPUT_CCA(1,1);
            TS1_full(I,4) = OUTPUT_CCA(1,5);
            TS1_full(I,5) = OUTPUT_CCA(1,2);
            TS1_full(I,6) = OUTPUT_CCA(1,3);
            TS1_full(I,7) = OUTPUT_CCA(1,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Fifth line for group netmats and amps, and subject spatial maps (GG-)
    A = fullfile(S,sprintf('results_%s_using_group_amps.mat',O{n}));
    if exist(A,'file');
        load(A)
        TS1_full(I,1) = OUTPUT(1,2);
        TS1_full(I,2) = OUTPUT(1,1);
        if exist('OUTPUT_CCA','var')
            TS1_full(I,3) = OUTPUT_CCA(1,1);
            TS1_full(I,4) = OUTPUT_CCA(1,5);
            TS1_full(I,5) = OUTPUT_CCA(1,2);
            TS1_full(I,6) = OUTPUT_CCA(1,3);
            TS1_full(I,7) = OUTPUT_CCA(1,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
end

%%%%% GET TABLE 3:
O = {'200','Yeo'};
T3_thresh = nan(12,7);
for n = 1:length(O);
    
    % thresholded
    if n == 1; I = 1; elseif n == 2; I = 7; end
    A = fullfile(S,sprintf('results_%s_using_group_amps_map_thresh_1.mat',O{n}));
    if exist(A,'file');
        load(A)
        T3_thresh(I,1) = OUTPUT(1,2); T3_thresh(I+3,1) = OUTPUT(2,2);
        T3_thresh(I,2) = OUTPUT(1,1); T3_thresh(I+3,2) = OUTPUT(2,1);      
        if exist('OUTPUT_CCA','var')
            T3_thresh(I,3) = OUTPUT_CCA(1,1); T3_thresh(I+3,3) = OUTPUT_CCA(3,1);
            T3_thresh(I,4) = OUTPUT_CCA(1,5); T3_thresh(I+3,4) = OUTPUT_CCA(3,5);
            T3_thresh(I,5) = OUTPUT_CCA(1,2); T3_thresh(I+3,5) = OUTPUT_CCA(3,2);
            T3_thresh(I,6) = OUTPUT_CCA(1,3); T3_thresh(I+3,6) = OUTPUT_CCA(3,3);
            T3_thresh(I,7) = OUTPUT_CCA(1,4); T3_thresh(I+3,7) = OUTPUT_CCA(3,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
          
    % binarised
    if n == 1; I = 2; elseif n == 2; I = 8; end
    A = fullfile(S,sprintf('results_%s_using_group_amps_map_thresh_1_bin.mat',O{n}));
    if exist(A,'file');
        load(A)
        T3_thresh(I,1) = OUTPUT(1,2); T3_thresh(I+3,1) = OUTPUT(2,2);
        T3_thresh(I,2) = OUTPUT(1,1); T3_thresh(I+3,2) = OUTPUT(2,1);      
        if exist('OUTPUT_CCA','var')
            T3_thresh(I,3) = OUTPUT_CCA(1,1); T3_thresh(I+3,3) = OUTPUT_CCA(3,1);
            T3_thresh(I,4) = OUTPUT_CCA(1,5); T3_thresh(I+3,4) = OUTPUT_CCA(3,5);
            T3_thresh(I,5) = OUTPUT_CCA(1,2); T3_thresh(I+3,5) = OUTPUT_CCA(3,2);
            T3_thresh(I,6) = OUTPUT_CCA(1,3); T3_thresh(I+3,6) = OUTPUT_CCA(3,3);
            T3_thresh(I,7) = OUTPUT_CCA(1,4); T3_thresh(I+3,7) = OUTPUT_CCA(3,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
    
    % Percentile binarised
    if n == 1; I = 3; elseif n == 2; I = 9; end
    A = fullfile(S,sprintf('results_%s_using_group_amps_map_thresh_-95_bin.mat',O{n}));
    if exist(A,'file');
        load(A)
        T3_thresh(I,1) = OUTPUT(1,2); T3_thresh(I+3,1) = OUTPUT(2,2);
        T3_thresh(I,2) = OUTPUT(1,1); T3_thresh(I+3,2) = OUTPUT(2,1);      
        if exist('OUTPUT_CCA','var')
            T3_thresh(I,3) = OUTPUT_CCA(1,1); T3_thresh(I+3,3) = OUTPUT_CCA(3,1);
            T3_thresh(I,4) = OUTPUT_CCA(1,5); T3_thresh(I+3,4) = OUTPUT_CCA(3,5);
            T3_thresh(I,5) = OUTPUT_CCA(1,2); T3_thresh(I+3,5) = OUTPUT_CCA(3,2);
            T3_thresh(I,6) = OUTPUT_CCA(1,3); T3_thresh(I+3,6) = OUTPUT_CCA(3,3);
            T3_thresh(I,7) = OUTPUT_CCA(1,4); T3_thresh(I+3,7) = OUTPUT_CCA(3,4);
        end
    end
    I = I+1;
    clear OUTPUT OUTPUT_CCA
end

%%%%% Table 1

O = {'025','200','PFM50','Glasser','Yeo'};
T1_CCA = nan(20,6);
for n = 1:length(O);
    I = (n-1)*4+1;
    % Real CCA results full netmats
    A = fullfile(S,sprintf('results_%s.mat',O{n}));
    if exist(A,'file');
        load(A)
        if exist('OUTPUT_CCA','var')
            T1_CCA(I,3) = OUTPUT_CCA(2,1);
            T1_CCA(I,4) = OUTPUT_CCA(2,5);
            T1_CCA(I,5) = OUTPUT_CCA(2,2);
            T1_CCA(I,6) = OUTPUT_CCA(2,6);
            T1_CCA(I,7) = OUTPUT_CCA(2,3);
        end
    end
    clear OUTPUT OUTPUT_CCA
    
    % Real CCA results partial netmats
    A = fullfile(S,sprintf('results_%s.mat',O{n}));
    if exist(A,'file');
        load(A)
        if exist('OUTPUT_CCA','var')
            T1_CCA(I+1,3) = OUTPUT_CCA(4,1);
            T1_CCA(I+1,4) = OUTPUT_CCA(4,5);
            T1_CCA(I+1,5) = OUTPUT_CCA(4,2);
            T1_CCA(I+1,6) = OUTPUT_CCA(4,6);
            T1_CCA(I+1,7) = OUTPUT_CCA(4,3);
        end
    end
    clear OUTPUT OUTPUT_CCA
end
    
% Real CCA results amplitudes
load('Results/results_CCA_table1_NEW.mat');
T1_CCA(3:4,3) = OUTPUT_CCA([5 6],1);
T1_CCA(3:4,4) = OUTPUT_CCA([5 6],2);
T1_CCA(3:4,5) = OUTPUT_CCA([5 6],3);
T1_CCA(3:4,6) = OUTPUT_CCA([5 6],4);
T1_CCA(3:4,7) = OUTPUT_CCA([5 6],5);

T1_CCA(7:8,3) = OUTPUT_CCA([3 4],1);
T1_CCA(7:8,4) = OUTPUT_CCA([3 4],2);
T1_CCA(7:8,5) = OUTPUT_CCA([3 4],3);
T1_CCA(7:8,6) = OUTPUT_CCA([3 4],4);
T1_CCA(7:8,7) = OUTPUT_CCA([3 4],5);

T1_CCA(9:12,3) = OUTPUT_CCA(7:10,1);
T1_CCA(9:12,4) = OUTPUT_CCA(7:10,2);
T1_CCA(9:12,5) = OUTPUT_CCA(7:10,3);
T1_CCA(9:12,6) = OUTPUT_CCA(7:10,4);
T1_CCA(9:12,7) = OUTPUT_CCA(7:10,5);

T1_CCA(15:16,3) = OUTPUT_CCA([1 2],1);
T1_CCA(15:16,4) = OUTPUT_CCA([1 2],2);
T1_CCA(15:16,5) = OUTPUT_CCA([1 2],3);
T1_CCA(15:16,6) = OUTPUT_CCA([1 2],4);
T1_CCA(15:16,7) = OUTPUT_CCA([1 2],5);

T1_CCA(19:20,3) = OUTPUT_CCA([12 11],1);
T1_CCA(19:20,4) = OUTPUT_CCA([12 11],2);
T1_CCA(19:20,5) = OUTPUT_CCA([12 11],3);
T1_CCA(19:20,6) = OUTPUT_CCA([12 11],4);
T1_CCA(19:20,7) = OUTPUT_CCA([12 11],5);
T1_CCA = T1_CCA(:,3:end);
clear OUTPUT OUTPUT_CCA

%%%%% GET NEW TABLE WITH REGRESSED RESULTS
O = {'025'};
Tnew_regressed = nan(8,7);
% Put in original results
Tnew_regressed(1,:) = TS1_full(10,:);
Tnew_regressed(3,:) = TS1_full(8,:);
Tnew_regressed(5,:) = T2_partial(10,:);
Tnew_regressed(7,:) = T2_partial(8,:);

% Driven by spatial only - netmats regressed out first
A = fullfile(S,sprintf('results_netsRegressed_%s_using_group_amps.mat',O{1}));
if exist(A,'file');
    load(A)
    Tnew_regressed(2,1) = OUTPUT(1,2); Tnew_regressed(6,1) = OUTPUT(2,2);
    Tnew_regressed(2,2) = OUTPUT(1,1); Tnew_regressed(6,2) = OUTPUT(2,1);
    if exist('OUTPUT_CCA','var')
        Tnew_regressed(2,3) = OUTPUT_CCA(1,1); Tnew_regressed(6,3) = OUTPUT_CCA(3,1);
        Tnew_regressed(2,4) = OUTPUT_CCA(1,5); Tnew_regressed(6,4) = OUTPUT_CCA(3,5);
        Tnew_regressed(2,5) = OUTPUT_CCA(1,2); Tnew_regressed(6,5) = OUTPUT_CCA(3,2);
        Tnew_regressed(2,6) = OUTPUT_CCA(1,3); Tnew_regressed(6,6) = OUTPUT_CCA(3,3);
        Tnew_regressed(2,7) = OUTPUT_CCA(1,4); Tnew_regressed(6,7) = OUTPUT_CCA(3,4);
    end
end
clear OUTPUT OUTPUT_CCA

% Driven by netmats only - spatial regressed out first
A = fullfile(S,sprintf('results_netsRegressed_%s_using_group_maps_using_group_amps_using_subject_netmats.mat',O{1}));
if exist(A,'file');
    load(A)
    Tnew_regressed(4,1) = OUTPUT(1,2); Tnew_regressed(8,1) = OUTPUT(2,2);
    Tnew_regressed(4,2) = OUTPUT(1,1); Tnew_regressed(8,2) = OUTPUT(2,1);
    if exist('OUTPUT_CCA','var')
        Tnew_regressed(4,3) = OUTPUT_CCA(1,1); Tnew_regressed(8,3) = OUTPUT_CCA(3,1);
        Tnew_regressed(4,4) = OUTPUT_CCA(1,5); Tnew_regressed(8,4) = OUTPUT_CCA(3,5);
        Tnew_regressed(4,5) = OUTPUT_CCA(1,2); Tnew_regressed(8,5) = OUTPUT_CCA(3,2);
        Tnew_regressed(4,6) = OUTPUT_CCA(1,3); Tnew_regressed(8,6) = OUTPUT_CCA(3,3);
        Tnew_regressed(4,7) = OUTPUT_CCA(1,4); Tnew_regressed(8,7) = OUTPUT_CCA(3,4);
    end
end
clear OUTPUT OUTPUT_CCA

save('Results/OUTPUT_TABLES.mat','T1_CCA','TS1_full','T2_partial','T3_thresh','Tnew_regressed');
    