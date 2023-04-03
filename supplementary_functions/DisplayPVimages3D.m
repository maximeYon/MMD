%% Compare two  immages
clearvars; close all; clc;

home_path = pwd; home_path = split(home_path,filesep);
home_path = join(home_path(1:end-2,1),filesep,1); home_path = home_path{1};
addpath(genpath([home_path filesep 'MATLABfunctions']));

%% Define path
pathData = 'C:\Users\Administrateur\Mon_Drive\Matlab\ProcessingPV360\data\20211220_tumor2\239\pdata\1';

% Import reconstructed images
imageObj = ImageDataObject(pathData);
img = squeeze(imageObj.data);
img = permute(img,[2 1 3 4]);
FoV = permute(imageObj.Visu.VisuCoreExtent,[2 1 3]);
x = 0:FoV(2)/(size(img,2)-1):FoV(2);
y = 0:FoV(1)/(size(img,1)-1):FoV(1);
size(img)
slice = 32;
cmaps = [0 max(img(:))/1.2];

% %% Display
% figure()
% set(gcf,'Position',[514.0000  217.0000  749.8966  574.0000])
% subplot(1,3,1)
% imagesc(y,x,squeeze(img(:,:,slice,1)),cmaps)
% set(gca,'Ydir','normal');
% axis image;
% title('Short TE')
% xlabel('FoV (mm)') ; ylabel('FoV (mm)')
% subplot(1,3,2)
% imagesc(y,x,squeeze(img(:,:,slice,2)),cmaps)
% set(gca,'Ydir','normal');
% axis image;
% title('b0 8 ms shapes')
% xlabel('FoV (mm)') ; ylabel('FoV (mm)')
% subplot(1,3,3)
% imagesc(y,x,squeeze(img(:,:,slice,3)),cmaps)
% set(gca,'Ydir','normal');
% axis image;
% title('b weighted 8 ms shapes')
% xlabel('FoV (mm)') ; ylabel('FoV (mm)')
% colormap gray

%% Create Gif
filename = 'RAREvfa_DOR166.gif';
h = figure();
set(gcf,'Position',[976.7241  229.4138  287.1725  561.5862])
set(gcf,'color','w');

for n = 1:size(img,4)
    
    imagesc(y,x,squeeze(img(:,:,slice,n)))
    set(gca,'Ydir','normal');
    axis image;
    title('RARE vfa DOR 166 protocol')
    xlabel('FoV (mm)') ; ylabel('FoV (mm)')
    colormap gray
    drawnow
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if n == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end








