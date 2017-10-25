clear all; close all; clc

% Load results
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
load('Results/CCAoutputs_withp.mat');
Ms = [CCAout(:,1) CCApart(:,1)];
Is = [CCAout(:,3) CCApart(:,3)];
names = strrep(names,'_',' ');

% Calculate confidence intervals
%Es1 = squeeze(max(grotRp(:,1,:))) - squeeze(prctile(grotRp(2:end,1,:),95));
%Es2 = squeeze(max(grotRpart(:,1,:))) - squeeze(prctile(grotRpart(2:end,1,:),95)); 
%Es = [Es1 Es2];
Es = [squeeze(std(grotRp(2:end,1,:))).*1.96 squeeze(std(grotRpart(2:end,1,:))).*1.96];

% Calculate average significance threshold for CCA results
AvSig = grotRp(2:end,1,:); AvSig = prctile(AvSig(:),95);
AvSigL = grotRp(2:end,1,[1:4 7:8]); AvSigL = prctile(AvSigL(:),95);

% Plot figure
N = {'Group-task-based rfMRI full netmat','Subject-task-based rfMRI full netmat'...
    'Group-task-based rfMRI partial netmat','Subject-task-based rfMRI partial netmat',...
    'Group-task-based rfMRI spatial maps','Subject-task-based rfMRI spatial maps'};

figure; subplot('position',[0.08 0.25 0.9 0.7]); set(gcf,'Position',[100 100 600 550],'PaperPositionMode','auto')
barweb(Ms,Es,1,names)
axis([0.5 6.5 0.5 0.8])
hold on; hline(AvSig,':k'); 
XYrotalabel([90 0])
legend({'CCA on complete data','Partial CCA'})
ylabel('r_U_V')

figure; subplot('position',[0.2 0.4 0.75 0.55]); set(gcf,'Position',[100 100 630 250],'PaperPositionMode','auto')
bar(0.5:2.5:6*2.5,CCAout([1:4 7:8],1),0.4,'b')
hold on
bar(1.5:2.5:6*2.5,CCApart([1:4 7:8],1),0.4,'r'); 
axis([-0.5 6*2.5+0.4 0.6 0.78]); A = hline(AvSigL,'-.k'); set(A,'LineWidth',2); vline(9.75,'k')
%ylabel('r');
set(gca,'xtick',1:2.5:6*2.5,'xticklabel',N); 
%set(gca,'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6 AvSigL 0.7 0.8],'yticklabel',{'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' sprintf('%1.2f',AvSigL) '0.7' '0.8'})
XYrotalabel([25 0])
legend({'full CCA','partial CCA'},'Location','NorthWest')
print(gcf,'-dtiff','-r300','Results/partial_CCA')




