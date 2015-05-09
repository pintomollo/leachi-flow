function install_leachi_flow
% INSTALL_LEACHI_FLOW adds the required directories to the matlabpath
% and handles the directories structure and the dependent libraries. It also
% compiles the required MEX libraries, hence a C/C++ compiler is requried.
%
% Wilson lab, University of Otago
% Simon Blanchoud
% 28.01.2015

  % Start by moving inside the leachi_flow folder
  cell_folder = which('install_leachi_flow');
  [current_dir, junk, junk] = fileparts(cell_folder);
  [root_dir, junk, junk] = fileparts(current_dir);
  cd(current_dir);

  % Add the proper directories to MATLAB path
  addpath(current_dir);
  addpath(fullfile(current_dir, 'GUI'));
  addpath(fullfile(current_dir, 'MEX'));
  addpath(fullfile(current_dir, 'file_io'));
  addpath(fullfile(current_dir, 'helpers'));
  addpath(fullfile(current_dir, 'image_analysis'));
  addpath(fullfile(current_dir, 'libraries'));
  addpath(fullfile(current_dir, 'pipeline'));
  addpath(fullfile(current_dir, 'simulation'));
  savepath;

  % And for the LOCI as well
  if (exist(fullfile(current_dir, 'bftools'), 'dir'))
    addpath(fullfile(current_dir, 'bftools'));
    savepath;
  end

  % Otherwise, try to insall it !
  if (exist('bfconvert.bat', 'file') ~= 2)
    button = questdlg('Should we try to install the Bio-Formats command line tools ?');

    % Ask for the user to confirm this foolness
    if (strncmpi(button, 'yes', 3))
      try
        rmdir('bftools', 's');
      catch
        % Nothing...
      end

      % This looks like a permanent link... up to now at least
      try
        disp('Downloading Bio-Formats tools...')

        unzip('http://downloads.openmicroscopy.org/latest/bio-formats/artifacts/bftools.zip');
        addpath(fullfile(current_dir, 'bftools'));
        savepath;
      catch
        errs = lasterror;
        warning('Tracking:installLOCI', ['Installation failed for the following reason:\n' errs.message]);
      end

      % Amazing enough !!
      if (exist('bfconvert.bat', 'file') == 2)
        disp('Done !');
        disp(' ');
      else
        disp('Failed... try to get the Bio-Formats command-line tools from http://www.openmicroscopy.org and place it in the "cast" folder');
        disp(' ');
      end
    end
  end

  % Try to insall it the FFMPEG library !
  if (~ispref('ffmpeg', 'exepath'))
    button = questdlg({'Should we try to locate the FFMPEG command line tools ?', ...
                      '(You need to have the library installed already)'});

    % Ask for the user to confirm this foolness
    if (strncmpi(button, 'yes', 3))
      ffmpegsetup();

      if (~ispref('ffmpeg', 'exepath'))
          disp('If you want to install the FFMPEG library, please visit www.ffmpeg.org');
      end
      disp(' ');
    end
  end

  % Check if the sparse 64bits flag is needed
  if ~isempty(strfind(computer(),'64'))
    mexopts = ' -largeArrayDims';
  else
    mexopts = '';
  end

  % Ask for the configuration only once
  did_setup = false;
  cd('MEX')

  % Try to compile the necessary MEX files
  if (exist('median_mex') ~= 3)
    try
      if (~did_setup)
        mex -setup;
      end
      eval(['mex' mexopts ' median_mex.c']);
      did_setup = true;
    catch ME
      cd(root_dir);
      error('Tracking:MEX', ['Could not compile the required MEX function!\n' ME.message]);
    end
  end
  if (exist('gaussian_mex') ~= 3)
    try
      if (~did_setup)
        mex -setup;
      end
      eval(['mex' mexopts ' gaussian_mex.c']);
      did_setup = true;
    catch ME
      cd(root_dir);
      error('Tracking:MEX', ['Could not compile the required MEX function!\n' ME.message]);
    end
  end
  if (exist('bilinear_mex') ~= 3)
    try
      if (~did_setup)
        mex -setup;
      end
      eval(['mex' mexopts ' bilinear_mex.c']);
      did_setup = true;
    catch ME
      cd(root_dir);
      error('Tracking:MEX', ['Could not compile the required MEX function!\n' ME.message]);
    end
  end
  if (exist('dsegment') ~= 3)
    try
      if (~did_setup)
        mex -setup;
      end
      eval(['mex' mexopts ' dsegment.cpp']);
      did_setup = true;
    catch ME
      cd(root_dir);
      error('Tracking:MEX', ['Could not compile the required MEX function!\n' ME.message]);
    end
  end
  if (exist('draw_gaussians_mex') ~= 3)
    try
      if (~did_setup)
        mex -setup;
      end
      eval(['mex' mexopts ' draw_gaussians_mex.c']);
      did_setup = true;
    catch ME
      cd(root_dir);
      error('Tracking:MEX', ['Could not compile the required MEX function!\n' ME.message]);
    end
  end
  cd(root_dir);

  % These folders are required as well
  if (~exist('TmpData', 'dir'))
    mkdir('TmpData');
  end
  if (~exist('export', 'dir'))
    mkdir('export');
  end

  % Confirm to the user that everything went fine
  disp('Installation successful !');

  % Gnu GPL notice
  fprintf(1, ['\nB. leachi blood flow plateform,  Copyright (C) 2015  Simon Blanchoud\n', ...
    'This program comes with ABSOLUTELY NO WARRANTY;\n', ...
    'This is free software, and you are welcome to redistribute it\n', ...
    'under certain conditions; read licence.txt for details.\n\n']);

  return;
end
