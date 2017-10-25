C = cov([X_true Y_true]);
C = C(1:size(X_true,2),size(X_true,2)+1:end);
[A,B,R0,U,V,stats] = canoncorr(X_true,Y_true); R0(1)

[u,s,v]=nets_svds((X_true')',0);
%s = diag(s); t = s > 1e-10; s(t) = sqrt(s(t)); s = diag(s);
%CX = u(:,t) * diag(sqrt(s(t)));
CX = u * s * u';
[u,s,v]=nets_svds((Y_true')',0);
CY = u * s * u';
%CX = cov(X_true')^0.5;
%CY = cov(Y_true')^0.5;

%%

poop=[];
for i=1:1000
XY = mvnrnd(zeros(size([X_true Y_true])),cov([X_true Y_true]));
X = XY(:,1:size(X_true,2)); Y = XY(:,size(X_true,2)+1:end);
X = CX*X; Y = CY*Y;
[A,B,R,U,V,stats] = canoncorr(X,Y);
poop=[poop R(1)]
end
