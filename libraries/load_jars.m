function load_jars()
% Based on Miji.m from Fiji, loads the jar files found in the "jars" forlder into memory.

  %% This script sets up the classpath to Fiji and optionally starts MIJ
  % Author: Jacques Pecreaux, Johannes Schindelin, Jean-Yves Tinevez

  %% Get the source directory (one up from this file)
  curr_directory = fileparts(fileparts(mfilename('fullpath')));

  %% Get the Java classpath
  classpath = javaclasspath('-all');

  %% Add all libraries in jars/ to the classpath

  % Switch off warning
  warning_state = warning('off');

  % Add the jars
  fprintf(1, ' Loading .jar files:     ')
  add_to_classpath(classpath, fullfile(curr_directory,'jars'));
  fprintf(1, '\b\b\b\bdone\n');

  % Switch warning back to initial settings
  warning(warning_state)

  return;
end

function add_to_classpath(classpath, directory)
  % Get all .jar files in the directory
  dirData = dir(directory);
  dirIndex = [dirData.isdir];
  jarlist = dir(fullfile(directory,'*.jar'));
  path_= cell(0);
  for i = 1:length(jarlist)
      fprintf('\b\b\b%3d', i);
      if not_yet_in_classpath(classpath, jarlist(i).name)
          path_{length(path_) + 1} = fullfile(directory,jarlist(i).name);
      end
  end

  %% Add them to the classpath
  if ~isempty(path_)
      javaaddpath(path_, '-end');
  end

  %# Recurse over subdirectories
  subDirs = {dirData(dirIndex).name};
  validIndex = ~ismember(subDirs,{'.','..'});

  for iDir = find(validIndex)
    nextDir = fullfile(directory,subDirs{iDir});
    add_to_classpath(classpath, nextDir);
  end

  return;
end

function test = not_yet_in_classpath(classpath, filename)
  %% Test whether the library was already imported
  expression = strcat([filesep filename '$']);
  test = isempty(cell2mat(regexp(classpath, expression)));

  return;
end
