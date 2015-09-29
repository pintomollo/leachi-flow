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

  for nc = 1:Nc

      filename = files{nc};
      [filepath, fname, fileext] = fileparts(filename);
      new_name = fullfile(out_path, [fname '_res' params.mesh_resolution fileext]);

      fprintf('Loading STL surface for channel %i...\n', nc);

      [face, vertex] = stlread(filename);

      fprintf('Simplifying the surface for channel %i...\n', nc);

      face = face.';
      vertex = vertex.';
      [face, vertex] = oocs_mex(face, vertex, params.mesh_resolution);
      face = face.';
      vertex = vertex.';

      fprintf('saving simplified surface for channel %i...\n', nc);

      stlwrite(fullfile(out_path, new_name, face, vertex);
  end
