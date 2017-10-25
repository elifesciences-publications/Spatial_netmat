%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code from steve: /vols/Data/HCP/Phase2/scripts900/DO_2b_groupPCA.m:

run /vols/Data/HCP/Phase2/scripts900/setup.m
load('/vols/Data/HCP/Phase2/group900/migp/migp_MSMALL.mat');
dPCA=4500;
MIGP.cdata=W(1:dPCA,:)';

%%%% MIGP-based dense connectome - with Wishart-rolloff to improve CNR and avoid Ring-of-Fire
 %i=4500; id=400; grot=[1:i]'; grot=(0.2+0.8*exp(-(grot/id).^2)).*(1-(0.5*grot/i)).^4;   figure; plot(grot);
EigS=sqrt(sum(MIGP.cdata.^2))';       % get sqrt(eigenvalues)
EigD=EigS.^2; EigD=EigD/max(EigD);    %	get normalised eigenvalues
EigDn1=round(length(EigD)*0.6);
EigDn2=round(length(EigD)*0.8);
g1=10000; g3=100000;                  % golden search for best-fitting spatial-DoF
for i=1:20
  g2=(g1+g3)/2;
  EigNull=iFeta([0:0.0001:5],round(1.05*length(EigD)),g2)';
  EigNull=EigNull*sum(EigD(EigDn2:EigDn2+10))/sum(EigNull(EigDn2:EigDn2+10));
  grot=EigD(1:EigDn1)-EigNull(1:EigDn1);
  if min(grot)<0,  g1=g2;  else,  g3=g2;  end;
  [i g1 g2 g3]
end
EigDc=EigD-EigNull(1:length(EigD));  % subtract null eigenspectrum
grot=smooth(EigDc(50:end),50,'loess'); i=min(find(abs(grot)<1e-5))+50-10;
EigDc(i:end)= (1-(1:(1+length(EigDc)-i))/(1+length(EigDc)-i)).^2 * grot(i-50);
  subplot(1,3,1);   hold off; plot(EigNull  ); hold on; plot(EigD,'g'); plot(EigDc,'r'); hold off;
EigSc=sqrt(abs(EigDc)); EigScn=(EigSc>0).*EigSc./(EigS/EigS(1));  	% get correction factors
  subplot(1,3,2); plot([EigScn  EigSc]); % plot([EigScn  EigSc grot400 grot400.*EigS/EigS(1) grot4000 grot4000.*EigS/EigS(1)]);
  subplot(1,3,3); plot([ EigS EigS.*EigScn]);

% create Wishart-adjusted MIGP spatial eigenvectors
MIGProw = MIGP.cdata * diag(EigScn);
  
% 2nd level varnorm of raw MIGP using sqrt(var(raw MIGP) - var(ROW MIGP))  then PCA
VN2=sqrt(var(MIGP.cdata,[],2) - var(MIGProw,[],2));  VN2=VN2/mean(VN2);
MIGPvn2 = MIGP.cdata ./ repmat(VN2,1,size(MIGProw,2));
[uu,ss,vv]=nets_svds(double(MIGPvn2),size(MIGPvn2,2)-1);

EigS=diag(ss);
EigD=EigS.^2; EigD=EigD/max(EigD);    %	get normalised eigenvalues
EigDn1=round(length(EigD)*0.6);
EigDn2=round(length(EigD)*0.8);
g1=10000; g3=100000;                  % golden search for best-fitting spatial-DoF
for i=1:20
  g2=(g1+g3)/2;
  EigNull=iFeta([0:0.0001:5],round(1.05*length(EigD)),g2)';
  EigNull=EigNull*sum(EigD(EigDn2:EigDn2+10))/sum(EigNull(EigDn2:EigDn2+10));
  grot=EigD(1:EigDn1)-EigNull(1:EigDn1);
  if min(grot)<0,  g1=g2;  else,  g3=g2;  end;
  [i g1 g2 g3]
end
EigDc=EigD-EigNull(1:length(EigD));  % subtract null eigenspectrum
grot=smooth(EigDc(50:end),50,'loess'); i=min(find(abs(grot)<1e-5))+50-10;
EigDc(i:end)= (1-(1:(1+length(EigDc)-i))/(1+length(EigDc)-i)).^2 * grot(i-50);
  subplot(1,3,1);   hold off; plot(EigNull  ); hold on; plot(EigD,'g'); plot(EigDc,'r'); hold off;
EigSc=sqrt(abs(EigDc)); EigScn=(EigSc>0).*EigSc./(EigS/EigS(1));  	% get correction factors
  subplot(1,3,2); plot([EigScn  EigSc ]);
  subplot(1,3,3); plot([ EigS EigS.*EigScn]);

% create Wishart-adjusted MIGP spatial eigenvectors
MIGProw = uu * ss * diag(EigScn);

clear ans EigD EigDc EigDn1 EigDn2 EigNull EigS EigSc EigScn g1 g2 g3 grot i MIGP MIGPvn2 ss uu VN2 vv W BATCH CIFTI SCRATCH WBC SUBJECTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set paths:
path(path,fullfile('/vols/Scratch/janineb','matlab','cifti-matlab'))
path(path,'/home/fs0/janineb/scratch/HCP/DMN/DMN_functions/')

% Load PFM group maps
PFMmaps = ft_read_cifti('/vols/Scratch/janineb/PROFUMO/PFMnew_S820_M50_Aug16_FinalModel.dtseries.nii');
PFMmaps = PFMmaps.dtseries; PFMmaps(isnan(PFMmaps(:,1))==1,:) = [];

% Create ground truth netmats using MIGProw and PFM group maps
PFMmaps = pinv(PFMmaps)';
GTnet = (PFMmaps'*MIGProw)*(MIGProw'*PFMmaps);
GTnet = (GTnet ./ repmat(sqrt(abs(diag(GTnet))),1,50)) ./ repmat(sqrt(abs(diag(GTnet)))',50,1); 
figure; imagesc(GTnet,[-0.5 0.5]); colorbar; title('''Ground truth'' netmat'); set(gcf,'Position',[20 20 850 700],'PaperPositionMode','auto')

% Save output
save('Ground_truth_netmat_PFMnew.mat','GTnet');

