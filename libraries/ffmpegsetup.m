function ffmpegsetup
%FFMPEGSETUP   Set up FFmpeg Toolbox
%   Run FFMPEGSETUP before using FFmpeg Toolbox for the first time. User
%   will be asked for the location of FFmpeg executable (binary) if it
%   cannot be located automatically.

% Copyright 2013 Takeshi Ikuma
% History:
% rev. - : (06-19-2013) original release
%        : (04-21-2015) Simon Blanchoud, modified the end to prevent empty
%                       preferences to be set

% dlgtitle = 'FFMPEG Toolbox Setup';
% ffmpegdir = fileparts(which(mfilename));

  % get existing config
  if ispref('ffmpeg','exepath')
     ffmpegexe = getpref('ffmpeg','exepath');
  else
     ffmpegexe = '';
  end

  % first check if it is arleady in the system path (should be the case for
  % both linux & mac
  ffmpegnames = 'ffmpeg';
  [fail,~] = system([ffmpegnames ' -version']);
  if ~fail
     ffmpegexe = ffmpegnames;
     setpref('ffmpeg','exepath',ffmpegexe); % save for later

      disp('FFMPEG   ...done.');
     return;
  end

  % if not found, ask user for the location of the file
  switch lower(computer)
     case {'pcwin' 'pcwin64'}
        filter = {'ffmpeg.exe','FFMPEG EXE file';'*.exe','All EXE files (*.exe)'};
     otherwise % linux/mac
        filter = {'*','All files'};
  end

  [filename, pathname] = uigetfile(filter,'Locate FFMPEG executable file',ffmpegexe);
  if numel(filename)>0 && filename(1)~=0 % cancelled
     ffmpegexe = fullfile(pathname,filename);
     if any(ffmpegexe==' '), ffmpegexe = ['"' ffmpegexe '"']; end

     % try
     [fail,msg] = system([ffmpegexe ' -version']);
     if fail || isempty(regexp(msg,'^ffmpeg','once'))
        ffmpegexe = '';
        warning('Tracking:ffmpegsetup','Invalid FFMPEG executable specified.');
     else
        setpref('ffmpeg','exepath',ffmpegexe);
        disp('FFMPEG   ...done.');
     end
  else
     warning('Tracking:ffmpegsetup','No valid FFMPEG executable found.');
  end

  return;
end
