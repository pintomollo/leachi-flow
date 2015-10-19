function simplify_mesh(params)

  %% parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if nargin==0
    params = get_struct('For3D');
  end

  if isempty(params.filename) || params.mesh_resolution==0, return, end % cancel by user

  files = params.filename;

  dir_out = '_3Dsmaller';

  [files, out_path] = get_filenames(files, dir_out);

  Nc = length(files);
  if (Nc == 0), disp('nada??'), return, end

  fprintf(' Computing simpler 3D volume :    ')

  for nc = 1:Nc

      fprintf('\b\b\b\b%3d,', nc);

      filename = files{nc};
      [filepath, fname, fileext] = fileparts(filename);
      new_name = fullfile(out_path, [fname '_res' num2str(params.mesh_resolution) fileext]);

      msg = sprintf(' loading surface...');
      fprintf(msg);
      %fprintf('Loading STL surface for channel %i...\n', nc);

      [face, vertex] = stlread(filename);

      fprintf([repmat('\b', 1, length(msg))]);
      msg = sprintf(' simplifying surface...');
      fprintf(msg);

      %fprintf('Simplifying the surface for channel %i...\n', nc);

      face = face.';
      vertex = vertex.';
      [face, vertex] = oocs_mex(face, vertex, params.mesh_resolution);
      face = face.';
      vertex = vertex.';

      fprintf([repmat('\b', 1, length(msg))]);
      msg = sprintf(' saving surface...');
      fprintf(msg);

      fprintf('saving simplified surface for channel %i...\n', nc);

      stlwrite(new_name, face, vertex);

      fprintf([repmat('\b', 1, length(msg))]);
  end

  fprintf('\b\b\b\b\bdone\n');

  return;
