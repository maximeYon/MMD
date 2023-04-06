function mdm_nii2pdf(nii_fn, pdf_fn, opt)
% function mdm_nii2pdf(nii_fn, pdf_fn, opt)
% 
% Converts one or several nii files in nii_fn to pdf files

if (nargin < 3), opt = []; end

% loop over cell array, if given
if (iscell(nii_fn))
    for c = 1:numel(nii_fn)
        [tmp_path, tmp_name] = msf_fileparts(nii_fn{c});
        if (~isempty(pdf_fn)), tmp_path = fileparts(pdf_fn); end
        tmp_pdf_fn = fullfile(tmp_path, [tmp_name '.pdf']);
        mdm_nii2pdf(nii_fn{c}, tmp_pdf_fn, opt);
    end
    return;
end

Iplot = mdm_nii_read(nii_fn);

figure(1)
clf;

axes('position',[0 0 1 1])

if size(Iplot,1) == 3
    Icol = zeros(size(Iplot,2),size(Iplot,3),3);
    for n = 1:3
        Icol(:,:,n) = squeeze(Iplot(n,:,:))';
    end
    imagesc(Icol/256);
else
    imagesc(Iplot');
    colormap('gray');
end

set(gca,'YDir','normal')
axis equal, axis tight, axis off
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
aspect = diff(ylim)/diff(xlim);

set(gcf, 'PaperPosition', 2*[0 0 [1 aspect]],'PaperSize', 2*[1 aspect]);
eval(['print ' pdf_fn ' -loose -dpdf'])

