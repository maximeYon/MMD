%% test opening EPI image

pathTestData = 'C:\Users\yon\Documents\Weizmann\PV6_SPEN\20180218_170228_SPEN_readout_offset_test_1_1\42\pdata\1\';
%
dirname = pathTestData(1:end-8);
if exist([dirname,'visu_pars'],'file')==2
    head0=textread([dirname,'visu_pars'],'%s');  
    VisuCoreDataOffs_size  = str2num(char(head0(strmatch('##$VisuCoreDataOffs=',head0')+1)));
    VisuCoreDataOffs  = str2num(char(head0(strmatch('##$VisuCoreDataOffs=',head0')+3:strmatch('##$VisuCoreDataOffs=',head0')+2+VisuCoreDataOffs_size)));
    VisuCoreDataSlope_size  = str2num(char(head0(strmatch('##$VisuCoreDataSlope=',head0')+1)));
    VisuCoreDataSlope  = str2num(char(head0(strmatch('##$VisuCoreDataSlope=',head0')+3:strmatch('##$VisuCoreDataSlope=',head0')+2+VisuCoreDataSlope_size)));
    VisuCoreWordType_line  = char(head0(strmatch('##$VisuCoreWordType=',head0')));
    VisuCoreWordType = VisuCoreWordType_line(1,21:end);
    VisuCoreByteOrder_line  = char(head0(strmatch('##$VisuCoreByteOrder=',head0')));
    VisuCoreByteOrder = VisuCoreByteOrder_line(1,22:end);
end





imageObj = ImageDataObject(pathTestData);

% image = imageObj.data;

path = 'C:\Users\yon\Documents\Weizmann\Matlab';


% writes an Image readable by ParaVision to disk: genertes a 2dseq and a visu_pars file. obj.writeImage

% [varargin, exportVisu] =bruker_addParamValue(varargin, 'overrideExportVisu', '@(x) isstruct(x)', []);
% if isempty(exportVisu)
%     exportVisu=obj.exportVisu;
% end
% 
% if ~exist(exportVisu.imagewrite_path,'file')
%     mkdir(exportVisu.imagewrite_path);
% end
[VisuCoreData]=bruker_write2dseq( [path, filesep, '2dseq'], imageObj.exportVisu, imageObj.data );
obj.exportVisu.VisuCoreDataMax=VisuCoreData.VisuCoreDataMax;
obj.exportVisu.VisuCoreDataMin=VisuCoreData.VisuCoreDataMin;
bruker_writeVisu( exportVisu.imagewrite_path, exportVisu);
