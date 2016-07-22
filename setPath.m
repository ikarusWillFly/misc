function setPath(cfg)
%---------------------------------------
% Sets the path for the default directories:
% - <HD>\projects\<user>\<project>
% - <HD>\projects\<user>\toolboxes
%
% Inputs and defaults:
%   - cfg.project   = 'testing';
%   - cfg.hd        = 'C:\';
%   - cfg.user      = 'voraco';
%   - cfg.sbj       = 'TEST';
%---------------------------------------
global PTH
%%% input checking
field = 'project';            value = 'testing';
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'hd';                 value = 'C:\';
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'user';               value = 'voraco';
if ~isfield(cfg,field), cfg.(field) = value; end
field = 'sbj';                value = 'TEST';
if ~isfield(cfg,field), cfg.(field).code = value; end

%%% general paths
PTH.main   = [cfg.hd,'projects\'];
PTH.user   = [PTH.main,cfg.user,'\'];
PTH.tbx    = [PTH.user,'toolboxes\'];
PTH.study  = [PTH.user,cfg.project,'\'];
PTH.code   = [PTH.study,'code\'];
PTH.data   = [PTH.study,'data\',cfg.sbj,'_',date,'\'];
PTH.util   = [PTH.user,'utilities\'];
PTH.sounds = [PTH.util,'sounds\'];

%%% toolboxes
PTH.pvt    = [PTH.tbx,'private_toolbox\'];
PTH.misc   = [PTH.pvt,'misc\'];
PTH.vh     = [PTH.pvt,'virtualHand\'];

PTH.ukt    = [PTH.tbx,'UKT_Toolbox_1_0_0\'];
PTH.ukt    = [PTH.tbx,'UKT_Toolbox_last\'];
PTH.ft_old = [PTH.tbx,'fieldtrip\fieldtrip-20121001\'];
PTH.ft_new = [PTH.tbx,'fieldtrip\fieldtrip-20130401\'];
PTH.bvr    = [PTH.tbx,'brainvision\'];
PTH.nc     = [PTH.tbx,'neuroConn\DC-STIMULATOR MC\RemoteControl\1.6.0'];
PTH.phi    = [PTH.tbx,'phidgets\'];
PTH.ptb    = [PTH.tbx,'Psychtoolbox\'];
PTH.r2p    = [PTH.tbx,'rehastim2plus\'];
PTH.sbox   = [PTH.tbx,'switchbox\'];
PTH.lsl    = [PTH.tbx,'labstreaminglayer\'];
end