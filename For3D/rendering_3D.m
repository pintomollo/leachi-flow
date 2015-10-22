function params = rendering_3D(params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function rendering_3D(params)
%
%% This function is the main one for Full Organ Rendering in 3D (For3D)
%% from a complete collection of immunohistologically labelled slices
%
%% Preliminary steps, in ImageJ
% 1/ adjust contrast & remove background for each color if necessary
% 2/ correct/remove artifacts if necessary (flipped sections, fat tissue, etc.)
% 3/ run the StackReg plugin (Thevenaz et al, EPFL)
% 4/ crop & save as multi-tif
%
% Color thresholds may be assessed in ImageJ from montage image
% (=> seing all images of stack at once)
%
%% -- input --
% params: parameters can be set either using the GUI or directly by the command line:
% i.e.
% rendering_3D({'*.tif'; '0'; '[130 -1 -1]'})
% => testing default threshold values for green and blue channels
%
% rendering_3D({'*.tif'; '0'; '[130 35 200]'})
% => using manually set thresholds
%
% default:
% rendering_3D({'*.tif';'0';'[-1 -1 -1]';'6.5';'22';'0.2';'[1/3 1/20 1/20]';'[0 0]';'0'})
%
%% -- output --
% For an input stack named 'input_file.tif', 3D images (as well as intermediate steps)
% are automatically saved in a subfolder named'temp/input_file_3D', with the extension
% '_3D_NN.png' (NN ranging from 0 to 45).
%
% see also filter_stack, equalize_stack, FOR3D_dialog_box
%
%% (c) Arnauld SERGE, 2011
% Contact: arnauld.serge@univ-mrs.fr
%
%% See ref:
%% J Immunol. 2013 Jan 15;190(2):586-96. doi: 10.4049/jimmunol.1200119. Epub 2012 Dec 17
%
%% Three-dimensional visualization of the mouse thymus organization in health and immunodeficiency.
%% Irla M, Guenot J, Sealy G, Reith W, Imhof BA, Sergé A.
%% Department of Pathology and Immunology, Faculty of Medicine, University of Geneva, Switzerland
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if nargin==0
    params = get_struct('For3D');
  end

  if isempty(params.filename), return, end % cancel by user

  files = params.filename;
  thresholds = params.thresholds;

  dir_out = '_3D';

  [files, out_path] = get_filenames(files, dir_out);

  Nc = length(files);
  if (Nc == 0), disp('nada??'), return, end

  new_names = files;

  %% thresholds auto
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(' Finding intensity thresholds :     ');
  for nc = 1:Nc
      fprintf('\b\b\b%3d', nc);

      [stk, junk, type] = load_sparse_stack(files{nc}, params.sparse_thresholds(nc));
      [Nx, Ny, Nz] = size(stk);

      if thresholds(nc) < 0
          scaling = abs(thresholds(nc));

          %fprintf('finding automatic threshold %g. Iteration...', nc)

          T1 = mygraythresh(stk, type);
          T2 = full(opthr(reshape(stk, Nx, [])));
          
          if Nc>1
              switch nc % tricky combination, might be optimized for a given dataset, or set manually
                  case 1 % red
                      if params.detect_IHC
                          a1 = 0.3; a2 =0.7;
                      else % intermediate threshold to get only medulla, avoiding non spe cortex
                          a1 = 0.7; a2 =0.3;
                      end
                      T = min((a1*T1 + a2*T2), 255);
                      %fprintf('Red threshold computation: T1 = %.1f, T2 = %.1f, %g*T1 + %g*T2 = %.1f\r', T1, T2, a1, a2, T)
                      thresholds(nc) = T;
                      
                  case 2 % green 
                      if params.detect_IHC
                          a1 = 0.3; a2 =0.7; a3 = 2;
                          T = min((a1*T1 + a2*T2)*a3, 255);
                      %    fprintf('Green threshold computation:  T1 = %.1f, T2 = %.1f, (%g*T1 + %g*T2)*%g = %.1f\r', T1, T2, a1, a2, a3, T)
                      else % low threshold to get all staining
                          T = T2; %/2;
                      %    fprintf('Green threshold computation:  T = %.1f\r', T)
                      end
                      thresholds(nc) = T;
                      
                  case 3 % blue 
                      if params.detect_IHC
                          T = T2;
                     %     fprintf('Blue threshold computation:  T = %.1f\r', T2)
                      else % high threshold to get only capsule
                          a1 = 0.8; a2 =0.2; a3 = 1.5;
                          T = min((a1*T1 + a2*T2)*a3, 255); % (0.8*T1 + 0.2*T2)*1.5;
                     %     fprintf('Blue threshold computation:  T1 = %.1f, T2 = %.1f, (%g*T1 + %g*T2)*%g = %.1f\r', T1, T2, a1, a2, a3, T)
                      end
                      thresholds(nc) = T;
              end
          else % gray, Nc = 1
              thresholds = mean([T1 T2]);
          end

          thresholds(nc) = scaling * thresholds(nc);
      end % if thresholds(nf, nc) == -1
  end % for nc = 1:Nc
  fprintf('\b\b\b\bdone\n');

  %% ****** generate 3D view ******
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf(' Computing 3D volume :      ')

  pixels = [params.pixel_size params.pixel_size(1)*ones(1,2-length(params.pixel_size)) params.slice_width];

  %%  plot isosurf 
  for nc = 1:Nc

      fprintf('\b\b\b\b%3d,', nc);

      msg = sprintf(' generating surface...');
      fprintf(msg);
      %fprintf('generating surface for channel %i...\n', nc)

      [stk, junk, type] = load_sparse_stack(files{nc}, params.sparse_thresholds(nc));

      [face, vertex] = MarchingCubes(stk, thresholds(nc));

      clear stk junk type;

      vertex = bsxfun(@times, vertex, pixels);
      vertex = bsxfun(@minus, vertex, mean(vertex));

      new_name = fullfile(out_path, ['volume_iso' num2str(round(thresholds(nc))) '_c' num2str(nc, '%02i') '.stl']);

      fprintf([repmat('\b', 1, length(msg)) repmat(' ', 1, length(msg)) repmat('\b', 1, length(msg))]);
      msg = sprintf(' saving surface...');
      fprintf(msg);

      %fprintf('saving surface for channel %i...\n', nc);

      stlwrite(new_name, face, vertex);
      new_names{nc} = new_name;

      fprintf([repmat('\b', 1, length(msg)) repmat(' ', 1, length(msg)) repmat('\b', 1, length(msg))]);
  end
  fprintf('\b\b\b\b\bdone\n');

  params.thresholds = thresholds;
  params.filename = new_names;

  return

