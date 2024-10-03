clear all

Nvector = 10:1000;
UDSRcell = cell(length(Nvector),1);
tic
parfor n = 1:length(Nvector);
%for n = 1:length(Nvector);
    N = Nvector(n);
    
    % Uniformly distribute 162 particles across the surface of the unit sphere 
    [V,Tri,~,Ue]=ParticleSampleSphere('N',N); % this operation takes ~8 sec on my machine (6GB RAM & 2.8 GHz processor)

%     % Visualize optimization progress 
%     figure, plot(0:numel(Ue)-1,Ue,'.-') 
%     set(get(gca,'Title'),'String','Optimization Progress','FontSize',20) 
%     xlabel('Iteration #','FontSize',15) 
%     ylabel('Reisz s-energy','FontSize',15)


    %%
    X = V(:,1);
    Y = V(:,2);
    Z = V(:,3);
    TRI = Tri;

%     figure(1), clf
%     h = trimesh(TRI,X,Y,Z);
%     set(h,'EdgeColor','b')
%     axis equal
%     title(['N = ' num2str(N)])
%     pause(.1)

%     UDSR.N = length(X);
%     UDSR.x = X;
%     UDSR.y = Y;
%     UDSR.z = Z;
%     UDSR.w = ones(UDSR.N,1);
%     UDSR.tri = TRI;
%     UDSR.theta = acos(UDSR.z);
%     UDSR.phi = atan2(UDSR.y,UDSR.x);

    UDSR = struct('N',length(X),'x',X,'y',Y,'z',Z,'w',1,'tri',TRI,'theta',acos(Z),...
        'phi',atan2(Y,X));
    UDSRcell{n,1} = UDSR;
end
toc
%%
for n = 1:length(Nvector);
    UDSR = UDSRcell{n,1};

    figure(1), clf
    h = trimesh(UDSR.tri,UDSR.x,UDSR.y,UDSR.z);
    set(h,'EdgeColor','b')
    axis equal
    title(['N = ' num2str(UDSR.N)])
    pause(.1)

    eval(['save UDSRTriN' num2str(UDSR.N) ' UDSR'])
end

