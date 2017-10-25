clear all; close all; clc

addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/

load('Results/results_CCA_table1_CIs.mat'); Aci = A1;
load('Results/results_CCA_table1_thresh05.mat');
load('Results/results_CCA_table1_NEWER.mat');
setdiff(A1,Aci)

A1 = strrep(A1,'_',' ');
A1 = strrep(A1,'amp','amplitude');
A1 = strrep(A1,'spatial','spatial maps');
A1 = strrep(A1,'Task','Task contrast spatial');
A1 = strrep(A1,'Fnetmat','full network matrix');
A1 = strrep(A1,'Pnetmat','partial network matrix');
A1 = strrep(A1,'area','surface area');
A1 = strrep(A1,'Glasser maps','Glasser spatial maps');
A1 = strrep(A1,'PFM441','PFM50');
A1 = strrep(A1,'amplitudelitude','amplitude');

I = [18;19;5;6;20;21;3;4;7;8;9;10;24;25;12;11];
ci = CIs(I,:); a = A1(I); t = thresh05(I,4); 
figure
set(gcf,'Position',[100 100 800 600],'PaperPositionMode','auto')
subplot('Position',[0.25 0.09 0.75 0.86])
c = [repmat([0 0 0.7],4,1);...
    repmat([0 0 1],4,1);...
    repmat([0 0.5 0.5],4,1);...
    repmat([0 0.7 0.7],3,1);...
    [0 1 1]];
barweb(ci(:,1), (ci(:,4)-ci(:,3))./2, 1,[],'CCA results',[],[],c)
ylabel('R_u_v');
hold on; X = 0.601:0.8/length(I):1.43;
for n = 1:length(I)
    N = nan(1,length(I)+1); N(n:n+1) = repmat(t(n),1,2);
    plot(X,N,'r','linewidth',2);
    %text(mean(X([n n+1])),OUTPUT_CCA(I(n),1)+0.03,sprintf('%1.2f',OUTPUT_CCA(I(n),3)));
end
axis([0.6 1.4 0 0.85])
view([90 90])
set(gca,'xtick',0.62:0.8/length(I):1.38,'xticklabel',a,'fontsize',10)
print(gcf,'-dtiff','-r300','Results/table1_fig.tif')

Uall = zeros(819,length(I)-1);
for n = 1:length(I)-1
    A = Us.(Aci{I(n)});
    Uall(:,n) = A(:,1);
end
figure; imagesc(abs(corr(Uall)),[-1 1]); colorbar
set(gca,'ytick',1:length(I)-1,'yticklabel',a(1:end-1),'fontsize',10);

I = [22;23;1;2;14;15;16;17;13];
ci = CIs(I,:); a = A1(I); t = thresh05(I,4);
figure
set(gcf,'Position',[900 100 800 300],'PaperPositionMode','auto')
subplot('Position',[0.25 0.15 0.75 0.77])
c = [repmat([0 0.7 0],4,1);....
    repmat([0 0.5 0.5],4,1);...
    repmat([0 0.9 0],1,1)];
barweb(ci(:,1), (ci(:,4)-ci(:,3))./2, 0.8,[],'CCA results (reduced subject sample n=441)',[],[],c)
ylabel('R_u_v');
hold on; X = 0.601:0.8/length(I):1.43;
for n = 1:length(I)
    N = nan(1,length(I)+1); N(n:n+1) = repmat(t(n),1,2);
    plot(X,N,'r','linewidth',2);
    %text(mean(X([n n+1])),OUTPUT_CCA(I(n),1)+0.03,sprintf('%1.2f',OUTPUT_CCA(I(n),3)));
end
axis([0.6 1.4 0 0.95])
view([90 90])
set(gca,'xtick',0.63:0.8/length(I):1.38,'xticklabel',a,'fontsize',10)

Uall = zeros(441,length(I)+1);
for n = 1:length(I)
    A = Us.(Aci{I(n)});
    Uall(:,n) = A(:,1);
end
Uall(:,n+1) = A(:,2);
figure; imagesc(abs(corr(Uall)),[-1 1]); colorbar
set(gca,'ytick',1:length(I)+1,'yticklabel',a,'fontsize',10);

