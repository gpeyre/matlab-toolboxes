function [NFV,smf_fname] = qslim(FV,varargin);
%QSLIM - Mesh simplification, wrapper function for Garland's QSLIM executable program
% function [NFV,smf_fname] = qslim(FV,varargin);
% varargin should be entered in pairs '<option>','<arg>'.
% Valid pairs used in this wrapper are:
%
%  '-t', <n>
%
%   Specify the desired number of faces in the simplified model.  You are not
%   guaranteed to get exactly the number you ask for, but it should generally
%   be within 1-2 faces of the target.
%
%  '-q',[]
%
%   Run quietly.  QSlim will not output the usual status information.
%
%  '-O', <n> (capital O, order)
%
%   Specify the policy for selecting the target position for edge
%   contractions.  The following levels of optimization are available:
%
%        3 -- Pick point which minimizes error [default].
%        2 -- Pick best point along edge.
%        1 -- Pick best of endpoints or midpoint.
%        0 -- Pick best of the two endpoints.
%
%  '-m', <penalty>
%
%     Set the penalty for bad meshes (default is 1)
%
%  '-c', <ratio>
%
%     Set the desired compactness ratio (default is 0)
%
%  '-B', <weight>
%
%     Specifies the weight assigned to boundary constraint planes.  The default
%     value is currently 1000.  Specify a weight of 0 to disable boundary
%     constraints entirely.
%
%  '-W', <n>
%
%     Select the quadric weighting policy.  Available policies are:
%
%        0 -- Weight all quadrics uniformly
%        1 -- Weight by area of contributing triangle [default]
%        2 -- Weight by angle around vertex
%
%  '-j',[]
%
%     Only perform contractions that do not remove any faces.  You can use this
%     feature in conjunction with custom edge sets (see below) to effect a form
%     of stitching.  This is not terribly reliable, it's just an experimental
%     feature.
%
%  '-F',[]
%
%     This will cause the simplification algorithm to use iterative face
%     contraction rather than iterative edge contraction.  Generally speaking,
%     this has the following effects:
%
%        - Simplification is faster
%        - Less memory is consumed
%        - Geometric quality of the results is reduced
%
%  '-o',<fname>
%
%     The output simplified file, default is temp_qslim_out.smf.
%
% Wrapper option:
%
%    If FV is a string, then that smf file is read, rather than writing out FV to
%    a file. Use the -save option below to write the structure FV to a file.
%    Saves execution time on subsequent passes
%
%  '-save', <fname>,
%
%     The filename to which to write the structure FV. Default is temp_qslim.smf.
%
% Defaults are set to -t, 2000, -m, 1000, -o temp_qslime_out.smf
%
% See http://graphics.cs.uiuc.edu/~garland/software/qslim.html
%  http://graphics.cs.uiuc.edu/~garland/research/thesis.html
%  for details on QSLIM and dissertation with technical details.
%
% Examples
%
% reduce to 2000 faces, with mesh penalty of 1000
% NFV = qslim(FV);
%
% write the FV to the brain.smf file, return 10,000 faces
% NFV = qslim(FV,'-save','brain.smf','-t',10000);
%
% load the brain.smf file for the FV information, quiet mode, compactness ratio
% of 1
% NFV = qslim('brain.smf','-q',[],'-c',1);
%
%
% See also REDUCEPATCH

%<autobegin> ---------------------- 27-Jun-2005 10:45:28 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - General
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\load_smf.m
%   toolbox\save_smf.m
%
% At Check-in: $Author: Mosher $  $Revision: 7 $  $Date: 6/27/05 9:00a $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005
%
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
%
% See BrainStorm website at http://neuroimage.usc.edu for further information.
%
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:45:28 -----------------------


% QSLIM program by
%
% Michael Garland
% Department of Computer Science
% University of Illinois
% 201 North Goodwin Avenue
% Urbana, IL 61801-2302
%
% Wrapper author:
% John C. Mosher

% ----------------------------- Script History ---------------------------------
% JCM 24-May-2004  Creation
% ----------------------------- Script History ---------------------------------

DEFAULT_INPUT = 'temp_qslim.smf'; % unless the user specifies in FV

if isstr(FV),
    % user specified an input file instead of a FV structure
    InputFile = FV;
else
    InputFile = DEFAULT_INPUT; % default
end

if ~ispc,
    error(sprintf('Unable to run QSLIM from this %s computer.',computer));
end

% where is the QSLIM on this computer
QSLIM = which('qslim.exe');

if isempty(QSLIM),
    error(sprintf('Cannot find qslim.exe in the public toolbox'));
end

OPERATION = varargin(1:2:end);
ARGUMENT = varargin(2:2:end);

if strmatch('-q',OPERATION),
    VERBOSE = 0; % silent running requested
else
    VERBOSE = 1;
end

% setup defaults
OutputStr = struct('option',[],'arg',[]);

OutputStr(1).option = '-t'; % wrapper argument
OutputStr(1).arg = 2000; % number of triangles generated

OutputStr(2).option = '-o'; % wrapper argument
OutputStr(2).arg = 'temp_qslim_out.smf'; % default

OutputStr(3).option = '-m';
OutputStr(3).arg = 1000; % helps prevent mesh badness

for i = 1:length(OPERATION),
    % have we already set this output?
    ndx = strmatch(OPERATION{i},{OutputStr.option});
    if isempty(ndx),
        % haven't defined
        % is it valid?
        switch OPERATION{i}
            case {'-O','-B','-W','-t','-F','-o','-I','-m','-c','-r','-M','-q','-j','-h'}
                % okay
                OutputStr(end + 1).option = OPERATION{i};
                OutputStr(end).arg = ARGUMENT{i};
            case '-save'
                InputFile = ARGUMENT{i}; % change the input file specification
            otherwise
                disp(sprintf('unknown qslim option %s',OPERATION{i}))
        end

    else
        % already defined, replace
        OutputStr(ndx).arg = ARGUMENT{i};
    end
end

% now make command string

commandstr = sprintf('"%s"',QSLIM); % initiate the start

for i = 1:length(OutputStr),
    switch OutputStr(i).option
        case {'-o','-I'}
            % string arg
            commandstr = sprintf(' %s %s %s',commandstr,OutputStr(i).option,OutputStr(i).arg);
        case {'-q','-j','-F','-r','-h'}
            % empty argument
            commandstr = sprintf(' %s %s',commandstr,OutputStr(i).option);
        case {'-B','-W','-t','-m'}
            % assume integer argument
            commandstr = sprintf(' %s %s %.0f',commandstr,OutputStr(i).option,OutputStr(i).arg);
        case {'-c'}
            % floating point argument
            commandstr = sprintf(' %s %s %f',commandstr,OutputStr(i).option,OutputStr(i).arg);
    end
end

% where was the output?
ndx = strmatch('-o',{OutputStr.option});
OutputFile = OutputStr(ndx).arg;

commandstr = sprintf('%s %s',commandstr,InputFile); % the input file

if ~isstr(FV),
    if VERBOSE
        disp(sprintf('Writing out FV to SMF file %s. ..',InputFile))
    end
    save_smf(InputFile,FV);
end

if VERBOSE
    disp('Executing . . .')
    disp(commandstr)
end
dos(commandstr);
drawnow

if VERBOSE
    disp(sprintf('Reading results from %s . . .',OutputFile))
end

NFV = load_smf(OutputFile);


function FV = load_smf(fname);
%LOAD_SMF - Load a simply written SMF file into the faces vertices structure
% function FV = load_smf(fname);
%
% Each line is one of
% # comment
% f <i> <j> <k>
% v <x> <y> <z>
% otherwise the line is ignored
%
% See also LOAD_SMF

%<autobegin> ---------------------- 27-Jun-2005 10:44:59 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Data Processing
%
% At Check-in: $Author: Mosher $  $Revision: 7 $  $Date: 6/27/05 9:00a $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005  
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:44:59 -----------------------


% Author: John C. Mosher, Ph.D.

% ----------------------------- Script History ---------------------------------
% JCM 24-May-2004  Creation
% ----------------------------- Script History ---------------------------------

fid = fopen(fname,'rt');

% preallocate
FV.faces = zeros(200000,3);
FV.vertices = zeros(100000,3);
iV = 0; % initialize vertices counter
iF = 0; % faces counter
while 1
   tline = fgetl(fid);
   if ~ischar(tline),
      break
   end
   
   [OP,ARG] = strtok(tline);
   
   switch OP
      case {'begin','end'}
         % ignore
      case 'v'
         iV = iV + 1;
         FV.vertices(iV,:) = sscanf(ARG,'%g %g %g')';
      case 'f'
         iF = iF + 1;
         FV.faces(iF,:) = sscanf(ARG,'%f %f %f')';
      case '#'
         % comment line
         disp(tline);
      otherwise
         % do nothing
   end
end

fclose(fid);

FV.vertices = FV.vertices(1:iV,:);
FV.faces = FV.faces(1:iF,:);


function save_smf(fname,FV);
%SAVE_SMF - Save out a file in a simple form of the SMF format
% function save_smf(fname,FV);
% Writes out one line per vertex or face:
% v <x> <y> <z>
% f <i> <j> <k>
%
% See also LOAD_SMF

%<autobegin> ---------------------- 27-Jun-2005 10:45:38 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Data Processing
%
% At Check-in: $Author: Mosher $  $Revision: 7 $  $Date: 6/27/05 9:00a $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005  
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:45:38 -----------------------


% Author: John C. Mosher, Ph.D.

% ----------------------------- Script History ---------------------------------
% JCM 24-May-2004  Creation
% ----------------------------- Script History ---------------------------------


fid = fopen(fname,'wt');

for i = 1:size(FV.vertices,1),
   fprintf(fid,'v %.3f %.3f %.3f\n',FV.vertices(i,:));
end

for i = 1:size(FV.faces,1),
   fprintf(fid,'f %.0f %.0f %.0f\n',FV.faces(i,:));
end

fclose(fid);