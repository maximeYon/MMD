function SQ = fSuperQuad(DT,N,gamma)

DT.cl = (DT.lambda1-DT.lambda2)/(DT.lambda1+DT.lambda2+DT.lambda3);
DT.cp = 2*(DT.lambda2-DT.lambda3)/(DT.lambda1+DT.lambda2+DT.lambda3);
DT.cs = 3*DT.lambda3/(DT.lambda1+DT.lambda2+DT.lambda3);

eval(['load UDSRTriN' num2str(N)])
SQ = UDSR;
TR = TriRep(SQ.tri, SQ.x, SQ.y, SQ.z);
Nsubdiv = 3;
TR=SubdivideSphericalMesh(TR,Nsubdiv);
SQ.tri = TR.Triangulation;
SQ.verts = TR.X;
SQ.x = TR.X(:,1);
SQ.y = TR.X(:,2);
SQ.z = TR.X(:,3);
SQ.N = numel(SQ.x);
SQ.theta = acos(SQ.z);
SQ.phi = atan2(SQ.y,SQ.x);


SQ.gamma = gamma;

if DT.cl >= DT.cp
    SQ.alpha = (1-DT.cp)^SQ.gamma;
    SQ.beta = (1-DT.cl)^SQ.gamma;
    q.x = DT.lambda1*sign(cos(SQ.theta)).*abs(cos(SQ.theta)).^SQ.beta;
    q.y = DT.lambda2*(-1)*sign(sin(SQ.phi)).*abs(sin(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
    q.z = DT.lambda3*sign(cos(SQ.phi)).*abs(cos(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
elseif DT.cl < DT.cp
    SQ.alpha = (1-DT.cl)^SQ.gamma;
    SQ.beta = (1-DT.cp)^SQ.gamma;
    q.x = DT.lambda1*sign(cos(SQ.phi)).*abs(cos(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
    q.y = DT.lambda2*sign(sin(SQ.phi)).*abs(sin(SQ.phi)).^SQ.alpha.*sign(sin(SQ.theta)).*abs(sin(SQ.theta)).^SQ.beta;
    q.z = DT.lambda3*sign(cos(SQ.theta)).*abs(cos(SQ.theta)).^SQ.beta;
end

SQ_PAS.x = q.x;
SQ_PAS.y = q.y;
SQ_PAS.z = q.z;

R.gamma = [
cos(DT.gamma) -sin(DT.gamma) 0
sin(DT.gamma) cos(DT.gamma) 0
0 0 1];

R.beta = [
cos(DT.beta) 0 sin(DT.beta)
0 1 0
-sin(DT.beta) 0 cos(DT.beta)];

R.alpha = [
cos(DT.alpha) -sin(DT.alpha) 0
sin(DT.alpha) cos(DT.alpha) 0
0 0 1];

R.mat = R.gamma*R.beta*R.alpha;

SQ_LF.x = R.mat(1,1)*SQ_PAS.x + R.mat(1,2)*SQ_PAS.y + R.mat(1,3)*SQ_PAS.z;
SQ_LF.y = R.mat(2,1)*SQ_PAS.x + R.mat(2,2)*SQ_PAS.y + R.mat(2,3)*SQ_PAS.z;
SQ_LF.z = R.mat(3,1)*SQ_PAS.x + R.mat(3,2)*SQ_PAS.y + R.mat(3,3)*SQ_PAS.z;

SQ.verts = [SQ_LF.x SQ_LF.y SQ_LF.z];
%SQ.verts = [SQ_PAS.x SQ_PAS.y SQ_PAS.z];
%SQ.verts = [SQ.x SQ.y SQ.z];

TR = triangulation(SQ.tri, SQ.verts);
SQ.norms = vertexNormal(TR,(1:SQ.N)');

RGBval = .7*[1 1 1];

SQ.c = zeros(SQ.N,3);
SQ.c(:,1) = RGBval(1)*ones(SQ.N,1);
SQ.c(:,2) = RGBval(2)*ones(SQ.N,1);
SQ.c(:,3) = RGBval(3)*ones(SQ.N,1);

