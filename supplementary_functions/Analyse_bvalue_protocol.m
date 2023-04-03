%% Analyse Test_bvalue protocol
clearvars;close all;
% Add function to path
home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-1,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% List parameters
% 500 Mhz ex vivo
% bvalues = [0.454	1.076	2.102	3.633	5.769	8.611	12.261	16.820	22.387	29.065]; %  3 shapes
bvalues = [0.235	0.558	1.091	1.886	2.995	4.470	6.365	8.732	11.622	15.089]; %  4 shapes

% 7T in vivo 
% bvalues = [0.015	0.035	0.069	0.120	0.191	0.286	0.407	0.559	0.744	0.966	1.228	1.534	1.887]; % 4 shapes
% bvalues = [0.028	0.068	0.134	0.232	0.369	0.551	0.785	1.077	1.433	1.861	2.366	2.956	3.635]; % 3 shapes
% bvalues = [0.048	0.114	0.223	0.386	0.613	0.915	1.303	1.788	2.379	3.089	3.928	4.906	6.034]; % 2 shapes

nbvalue = size(bvalues,2);
Nfreq_shapes = 3;
Nshapes = 3;
Nb0 = 1;
N_noDiff = 1;

%% Path for data
data_path = 'optimization_exvivo_fresh\16'; 
% data_path = 'invivoRat\24'; 
procno =1;

data_path = split(data_path,'\'); data_path = join(data_path,filesep,1); data_path = data_path{1};
data_path = [home_path filesep 'data' filesep data_path];

%% Import data and list
% Import reconstructed images
imageObj = ImageDataObject([data_path filesep 'pdata' filesep '1']);
img = squeeze(imageObj.data);
Nsample = size(img,3);


%% Compare ROI data
figure(1)
set(gcf,'Position', [4 401 1863 577])
subplot(1,3,1)
imagesc(sum(abs(img),3));
colormap gray;
h = imfreehand;
mask = createMask(h);

LinData = img(repmat(mask,1,1,size(img,3)));
LinData = reshape(LinData,sum(mask(:)),size(img,3));
LinData = sum(LinData,1);
LinData = LinData./max(LinData(:)).*100;

%% Separate data per shape
b0 = 1+N_noDiff:(Nfreq_shapes*Nshapes+Nb0):Nsample;
for b = 1:nbvalue
    for sh = 1:Nshapes
        shape(sh,1+(b-1)*Nfreq_shapes:b*Nfreq_shapes)=1+N_noDiff+Nb0+(sh-1)*Nfreq_shapes+(b-1)*(Nfreq_shapes*Nshapes+Nb0):N_noDiff+Nb0+(sh)*Nfreq_shapes+(b-1)*(Nfreq_shapes*Nshapes+Nb0);
    end
end

for sh = 1:Nshapes
    ShapeLinData(sh,:) = LinData(1,shape(sh,:));
end

b0LinData = LinData(1,b0);
Xb0 = 1+N_noDiff:(Nfreq_shapes*Nshapes+Nb0):(Nfreq_shapes*Nshapes+Nb0)*nbvalue;
% Xshape = 3:5; 
Xshape = 1+N_noDiff+Nb0:N_noDiff+Nb0+Nfreq_shapes; 
for b = 1:nbvalue-1
    Xshape = [Xshape (1+N_noDiff+Nb0:N_noDiff+Nb0+Nfreq_shapes)+(Nfreq_shapes*Nshapes+Nb0)*b];
end

subplot(1,3,2:3)
hold on
plot(Xb0,b0LinData,'*k')
for sh = 1:Nshapes
    plot(Xshape,ShapeLinData(sh,:),'*')
end
legend('b0','sphere','plan','stick')
ylim([0 max(b0LinData)]);
xticks(Xb0);
xticklabels(string(bvalues))

