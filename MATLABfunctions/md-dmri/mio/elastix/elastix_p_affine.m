function p = elastix_p_affine(n_iter)
% function p = elastix_p_affine(n_iter)

if (nargin < 1), n_iter = 300; end

% ImageTypes
p.FixedInternalImagePixelType           = 'float';
p.FixedImageDimension                   = 3;
p.MovingInternalImagePixelType          = 'float';
p.MovingImageDimension                  = 3;

% Components setup
p.Registration                          = 'MultiResolutionRegistration';
p.FixedImagePyramid                     = 'FixedRecursiveImagePyramid';
p.MovingImagePyramid                    = 'MovingRecursiveImagePyramid';
p.Interpolator                          = 'BSplineInterpolator';
p.Optimizer                             = 'AdaptiveStochasticGradientDescent';
p.ResampleInterpolator                  = 'FinalBSplineInterpolator';
p.Resampler                             = 'DefaultResampler';
p.Transform                             = 'AffineDTITransform';
p.Metric                                = 'AdvancedMattesMutualInformation';

% Et c
p.ErodeMask                             = 'false';
p.HowToCombineTransforms                = 'Compose';
p.AutomaticTransformInitialization      = 'false';
p.AutomaticScalesEstimation             = 'true';
p.DefaultPixelValue                     = 0;

% Output in nifti
p.WriteResultImage                      = 'true';
p.ResultImagePixelType                  = 'float';
p.ResultImageFormat                     = 'nii';
p.CompressResultImage                   = 'false';

% Resolution and sampling
p.NumberOfResolutions                   = 1;
p.MaximumNumberOfIterations             = n_iter; 
p.NumberOfHistogramBins                 = [32];
p.ImagePyramicSchedule                  = [1 1 1];
p.ImageSampler                          = 'RandomCoordinate';
p.NewSamplesEveryIteration              = 'true';
p.NumberOfSpatialSamples                = 8192;

% Order of B-Spline interpolation used in each resolution level, final et c
p.BSplineInterpolationOrder             = 1;
p.FinalBSplineInterpolationOrder        = 3;

        

