%% Convert Python bootstrap into Matlab
close all; clearvars; clc;

%% Select bootstrap path
bs_path = 'C:\Users\User\Downloads\bootstraps\distributions\bootstraps';
bs_list = my_getdirno(bs_path);

%% Add function to path
addpath('npy-matlab-master\npy-matlab')

%% Unzip NPZ file
% unzip([bs_path filesep num2str(bs_list(1)) filesep 'mfs.npz'],[bs_path filesep num2str(bs_list(1)) filesep]);

%% Open NPY file
% py.importlib.import_module('numpy')
% test = py.numpy.load([bs_path filesep num2str(bs_list(1)) filesep 'fit_list.npy'],"True");
b = readNPY([bs_path filesep num2str(bs_list(1)) filesep 'fit_list.npy']);
% X = npyread([bs_path filesep num2str(bs_list(1)) filesep 'fit_list.npy']);

%% functions
function dirno = my_getdirno(dir_path)
dirlist = dir(dir_path);
[Ndir,~] = size(dirlist);
dirno = [];
for ndir = 1:Ndir
    dirnotest = str2num(getfield(dirlist,{ndir,1},'name'));
    if ~isempty(dirnotest)
        dirno = [dirno dirnotest];
    end
end
dirno = sort(dirno,'ascend');
end