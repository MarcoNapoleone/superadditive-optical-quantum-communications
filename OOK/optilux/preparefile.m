function [outfile,fid] = preparefile(dirname,file_name,issave)

%PREPAREFILE prepare output file for writing information
%   [OUTFILE,FID] = PREPAREFILE(DIRNAME,FILE_NAME,ISSAVE) creates a unique
%   output file OUTFILE if ISSAVE is true and returns also the FID of the
%   file identifier (see FOPEN). If ISASAVE=false, FID=1 and OUTFILE are 
%   set to the empty string.
%
%   This fuction is useful to simplify both the saving to file and printing
%   on the command window.
%
%   Author: Paolo Serena, 2021
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2021  Paolo Serena, <serena@tlc.unipr.it>
%
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if issave % save results to file?    
    if ~exist(dirname,'dir');
        mkdir(dirname);
    end
    outfile = [pwd,filesep,dirname,filesep,file_name]; % output file for savings
    nmot = 1; % if output file already exists, append counter
    while exist([outfile,'.txt'],'file')
        outfile = [pwd,filesep,dirname,filesep,file_name,'_',num2str(nmot)];
        nmot=nmot+1; % increase counter on output file name
    end    
    fid=fopen([outfile,'.txt'],'a');
else
    fid = 1; % do not save or write info to file
    outfile = '';
end
