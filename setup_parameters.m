% Setup_parameters for ambient noise processing
%
% NJA, 4/2/2016

%%% --- Paths to important files --- %%%
parameters.workingdir = pwd;
parameters.workingdir = [parameters.workingdir,'/'];

parameters.datapath = '~/Unix/Data/ambnoise/prep_data/';
% parameters.datapath = '/Users/accardo/Unix/Data/OBS/prep_data/';
% parameters.datapath = '/Volumes/Fraser/SEGMeNT_ambnoise/complete_array/prep_data/';
parameters.ccfpath = [parameters.workingdir,'ccf/'];
parameters.stalist = textread([parameters.workingdir,'stalist'],'%s\n');
parameters.nsta = length(parameters.stalist);

%%% --- Parameters for initial processing --- %%%
parameters.dt = 1; % sample rate
parameters.comp = 'HH'; % component
parameters.mindist = 20; % min. distance in kilometers

%%% --- Parameters for ccf_ambnoise --- %%%
parameters.winlength = 3;

%%% --- Parameters for fitbessel --- %%%
parameters.npts = parameters.winlength*3600;
