clear all; close all; clc

% Set up analysis
Nperm = 100000; 
addpath /vols/Scratch/janineb/HCP/CCA/files/

% Load inputs
load('Results/CCAinputs_taskStuff.mat');
names = fieldnames(CCAinput);
CCAout = zeros(length(names),4);
U = zeros(790,100,length(names)); V = zeros(790,100,length(names)); p = zeros(100,length(names)); grotRp = zeros(Nperm,101,length(names)); if Nperm>1; R = zeros(101,length(names)); else R = zeros(100,length(names)); end
CCApart = zeros(length(names),4);
Upart = zeros(790,100,length(names)); Vpart = zeros(790,100,length(names)); ppart = zeros(100,length(names)); grotRpart = zeros(Nperm,101,length(names)); if Nperm>1; Rpart = zeros(101,length(names)); else Rpart = zeros(100,length(names)); end

% Run CCA
for n = 1:length(names)
    A = CCAinput.(names{n});
    if Nperm>0 ;
        [~,U(:,:,n),V(:,:,n),R(:,n),Rmax,Imax,p(:,n),Pmax,grotRp(:,:,n)] = S_CCA(A,Nperm);
    else
        [~,U(:,:,n),V(:,:,n),R(:,n),Rmax,Imax] = S_CCA(A,Nperm);
    end
    CCAout(n,1) = R(1,n); CCAout(n,3) = Rmax;  CCAout(n,4) = Imax;
    if Nperm>1; CCAout(n,2) = p(1,n); end
    save('Results/CCAoutputs_withp.mat','CCAout','U','V','R','p','grotRp','CCApart','Upart','Vpart','Rpart','ppart','grotRpart','names');
end

% Partial CCA's
for n = 1:length(names)
    A = CCAinput.(names{n});
    if strfind(names{n},'grpmaps')
        B = CCAinput.Pgrp;
    elseif strfind(names{n},'submaps')
        B = CCAinput.Psub;
    elseif strfind(names{n},'grp')
        B = CCAinput.DR_groupmaps;
    elseif strfind(names{n},'sub')
        B = CCAinput.DR_subjectmaps;
    end
    A = A - B * pinv(B)*A;
    if Nperm>0 ;
        [~,Upart(:,:,n),Vpart(:,:,n),Rpart(:,n),Rmax,Imax,ppart(:,n),Pmax,grotRpart(:,:,n)] = S_CCA(A,Nperm);
    else
        [~,Upart(:,:,n),Vpart(:,:,n),Rpart(:,n),Rmax,Imax] = S_CCA(A,Nperm);
    end
    CCApart(n,1) = Rpart(1,n); CCApart(n,3) = Rmax;  CCApart(n,4) = Imax;
    if Nperm>1; CCApart(n,2) = ppart(1,n); end
    save('Results/CCAoutputs_withp.mat','CCAout','U','V','R','p','grotRp','CCApart','Upart','Vpart','Rpart','ppart','grotRpart','names');
end

% Test for significance differences between results
test = grotRp(:,1,1) - grotRp(:,1,2);
testp = (1+sum(test(2:end)>=test(1)))/Nperm;
fprintf('%1.5f\n',testp)
% test = grotRp(:,1,3) - grotRp(:,1,4);
% testp = (1+sum(test(2:end)>=test(1)))/Nperm;

