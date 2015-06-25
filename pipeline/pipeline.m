function pipeline(varargin)

  [myrecording, mysimulation, opts] = parse_inputs(varargin{:});

  return;
end

function [myrecording, mysimulation, opts] = parse_inputs(varargin)

  myrecording = get_struct('myrecording');
  mysimulation = get_struct('simulation');
  opts = get_struct('options');

  conf_rec = '';
  conf_sim = '';
  conf_opt = '';

  %%%% First run: get structures, load config files, remove varargins
  %%%% Second pass, load parameter pairs !

  skip = false;
  for i=1:length(varargin)
    if (skip)
      skip = false;
    else
      if (isstruct(varargin{i}))
        fitting = varargin{i};
      elseif (isnumeric(varargin{i}))
        if (isempty(fitting.init_pos))
          fitting.init_pos = varargin{i};
        else
          fitting.bounds = varargin{i};
        end
      elseif (ischar(varargin{i}))
        if (i == length(varargin))
          error('Properties pairs must come in PAIRS.');
        else
          % If the parameter exists in opts we simply assign it the
          % provided value
          if (isfield(fitting, varargin{i}))
            fitting.(varargin{i}) = varargin{i+1};

          % Or ignore it !
          else
            warning(['Property ''' varargin{i} ''' does not exist. Ignoring']);
          end

          skip = true;
        end
      elseif (isa(varargin{i}, 'function_handle'))
        fitting.error_function = varargin{i};
      else
        warning(['Unknown type of input ''' class(varargin{i}) '''. Ignoring']);
      end
    end
  end


  % Check what we got as inputs
  if (nargin > 0)
    mymovies = varargin{1};
    varargin(1) = [];

    if (length(varargin) > 0 & isstruct(varargin{1}))
      if (isfield(varargin{1}, 'reaction_params'))
        opts = varargin{1};
        varargin(1) = [];
      elseif (isfield(varargin{1}, 'parameter_set'))
        fitting = varargin{1};
        varargin(1) = [];
      end
      if (length(varargin) > 0 & isstruct(varargin{1}))
        if (isfield(varargin{1}, 'reaction_params'))
          opts = varargin{1};
          varargin(1) = [];
        elseif (isfield(varargin{1}, 'parameter_set'))
          fitting = varargin{1};
          varargin(1) = [];
        end
      end
    end

    % Now we check that the parameters were provided in pairs
    npairs = length(varargin) / 2;
    if (npairs ~= floor(npairs))
      uuid = str2double(varargin(1));
      varargin(1) = [];
      npairs = floor(npairs);
    end

    % Loop over the pairs of parameters
    for i = 1:npairs
      switch varargin{2*i - 1}
        case 'config_modeling'
          conf_mod = varargin{2*i};
        case 'config_fitting'
          conf_fit = varargin{2*i};
      end
    end

    if (~isempty(conf_mod))
      opts = load_parameters(opts, conf_mod);
    end

    if (~isempty(conf_fit))
      fitting = load_parameters(fitting, conf_fit);
    end

    % Loop over the pairs of parameters
    for i = 1:npairs

      % If the parameter exists in opts we simply assign it the
      % provided value
      if (isfield(fitting, varargin{2*i - 1}))
        fitting.(varargin{2*i - 1}) = varargin{2*i};
      elseif (isfield(opts, varargin{2*i - 1}))
        opts.(varargin{2*i - 1}) = varargin{2*i};
      else
        switch varargin{2*i - 1}
          case 'config_modeling'
            conf_mod = varargin{2*i};
          case 'config_fitting'
            conf_fit = varargin{2*i};
          otherwise
            warning on
            warning(['Property ''' varargin{2*i -1} ''' does not exist. Ignoring']);
        end
      end
    end
  end

  if (~fitting.fit_full)
    opts = load_parameters(opts, 'maintenance.txt');
  end

  if (~iscell(mymovies))
    mymovies = {mymovies};
  end

  return;
end
