function opt = qdti_opt(opt)
% function opt = qdti_opt(opt)

opt.qdti.present = 1;

opt.qdti = ensure_field(opt.qdti, 'tmp', 1); 