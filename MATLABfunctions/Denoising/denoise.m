function [denoisedImage,S2,P] = denoise(image,windowDim,mask,verboseFlag)
% The function takes as input arguments "image" containing the image data
% in an array with 3 or 4 dimensions. The last dimension must carry the
% variation with the magentic field, so that the first dimensions only
% correspond to differet pixels. "windowDim" specifies the sidelength of
% the pixel window defining the local neighborhood in which the denoising
% is carried out in each iteration. Depending on "image" having 3 or 4
% dimensions, "windowDim" should be a rowvector with 2 or 3 elements.
% "mask" is a logical array with dimensions corresponding to the first
% dimensions of "image", so every pixel is marked either for denoising
% (true) or skipping (false). "verboseFlag" is a boolean which when true
% makes the function output progress information to the matlab interface.
% Both "mask" and "verboseFlag" are optional inputs with the defaults
% verboseFlag=true and the mask set for denoising of the entire input
% image.
%
% This matlab function was implemented by Jonas Olesen and is based on the
% denoising algorithm for diffusion MRI data presented by Veraart et al.
% (2016).

%% adjust image dimensions and assert
if nargin<3 || numel(mask)==0
    imageDim = size(image);
    mask = true(imageDim(1:end-1));
end
if nargin<4
    verboseFlag = true;
end

imageDim = size(image);
imageDimOld = imageDim;
maskDim = size(mask);
if length(maskDim)==2, maskDim(3)=1; end
if length(imageDim)==3
    assert(all(maskDim(1:2)==imageDim(1:2)),'mask dimensions does not match image dimensions.')
    image = reshape(image,[imageDim(1:2) 1 imageDim(3)]);
    imageDim = [imageDim(1:2) 1 imageDim(3)];
else
    assert(length(imageDim)>2 && length(imageDim)<5,'incorrect dimensions of image data: must be of 3 or 4 dimensions.')
    assert(all(maskDim==imageDim(1:3)),'mask dimensions does not match image dimensions.')
end

if length(windowDim)==2
    windowDim(3) = 1;
else
    assert(length(windowDim)>1 && length(windowDim)<4,'window has incorrect dimensions.')
end
assert(all(windowDim<=imageDim(1:3)),'window exceeds image dimensions.')


%% denoise image
M = size(image,4);
N = prod(windowDim);
denoisedImage = zeros(size(image));
P = zeros(imageDim(1:3));
S2 = zeros(imageDim(1:3));
counter = zeros(imageDim(1:3));
progress0 = 0;
total = imageDim(1)-windowDim(1);
for i = 0:imageDim(1)-windowDim(1)
    for j = 0:imageDim(2)-windowDim(2)
        for k = 0:imageDim(3)-windowDim(3)
            % Check mask
            rows=i+(1:windowDim(1));
            cols = j+(1:windowDim(2));
            slis = k+(1:windowDim(3));
            maskCheck = reshape(mask(rows,cols,slis),[N 1])';
            if all(~maskCheck), continue, end
            
            % Create X data matrix
            X = reshape(image(rows,cols,slis,:),[N M])';
            
            % Remove voxels not contained in mask
            X(:,~maskCheck) = [];
            if size(X,2)==1, continue, end % skip if only one voxel of window in mask
            
            % Perform denoising
            newX=zeros(M,N); sigma2=zeros(1,N); p=zeros(1,N);
            [newX(:,maskCheck),sigma2(maskCheck),p(maskCheck)] = denoiseMatrix(X);
            if all(p==0), continue, end % skip if no signal carrying components detected
            
            % Assign newX to correct indices in denoisedImage
            denoisedImage(rows,cols,slis,:) = denoisedImage(rows,cols,slis,:) + reshape(newX',[windowDim M]);
            P(rows,cols,slis) = P(rows,cols,slis) + reshape(p,windowDim);
            S2(rows,cols,slis) = S2(rows,cols,slis) + reshape(sigma2,windowDim);
            counter(rows,cols,slis) = counter(rows,cols,slis)+1;
        end
    end
    if verboseFlag
        progress = round(i/total*10);
        if progress>progress0
            progress0 = progress;
            fprintf('row progress: %g%%\n',progress*10)
        end
    end
end
skipCheck = mask & counter==0;
counter(counter==0) = 1;
denoisedImage = bsxfun(@rdivide,denoisedImage,counter);
P = bsxfun(@rdivide,P,counter);
S2 = bsxfun(@rdivide,S2,counter);

%% adjust output to match input dimensions
% Assign original data to denoisedImage outside of mask and at skipped voxels
original = bsxfun(@times,image,~mask);
denoisedImage = denoisedImage + original;
original = bsxfun(@times,image,skipCheck);
denoisedImage = denoisedImage + original;

% Shape denoisedImage as orginal image
if length(imageDimOld)==3
    denoisedImage = reshape(denoisedImage,imageDimOld);
    S2 = reshape(S2,imageDimOld(1:end-1));
    P = reshape(P,imageDimOld(1:end-1));
end