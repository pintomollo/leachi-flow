function pipeline(varargin)

  [myrecording, mysimulation, opts] = parse_inputs(varargin{:});

  keyboard

  return;
end

function [myrecording, mysimulation, opts] = parse_inputs(varargin)

  myrecording = get_struct('myrecording');
  mysimulation = get_struct('simulation');
  opts = get_struct('options');

  conf_rec = '';
  conf_sim = '';
  conf_opt = '';

  % First run: get structures, load config files & remove them
  for i=length(varargin):-1:1
    if (isstruct(varargin{i}))
      if (isfield(varargin{i}, 'experiment'))
        myrecording = varargin{i};
        varargin(i) = [];
      elseif (isfield(varargin{i}, 'init_simulation'))
        mysimulation = varargin{i};
        varargin(i) = [];
      elseif (isfield(varargin{i}, 'config_files'))
        opts = varargin{i};
        varargin(i) = [];
      end
    elseif (ischar(varargin{i}))
      switch varargin{i}
        case 'config_recording'
          conf_rec = varargin{i+1};
          varargin(i:i+1) = [];
        case 'config_simulation'
          conf_sim = varargin{i+1};
          varargin(i:i+1) = [];
        case 'config_options'
          conf_opt = varargin{i+1};
          varargin(i:i+1) = [];
      end
    end
  end

  % Retrieve the latest version of the structures
  myrecording = update_structure(myrecording, 'myrecording');
  mysimulation = update_structure(mysimulation, 'simulation');
  opts = update_structure(opts, 'options');

  % Load the configuration files
  if (~isempty(conf_rec))
    myrecording = load_parameters(myrecording, conf_rec);
  end
  if (~isempty(conf_sim))
    mysimulation = load_parameters(mysimulation, conf_sim);
  end
  if (~isempty(conf_opt))
    opts = load_parameters(opts, conf_opt);
  end

  % Now we check that the parameters were provided in pairs
  npairs = length(varargin) / 2;
  if (npairs ~= floor(npairs))
    error('Properties pairs must come in PAIRS.');
  end

  % Second run: loop over the pairs of parameters and assign them
  for i = 1:npairs

    % Get the pair
    field = varargin{2*i - 1};
    value = varargin{2*i};

    % Try to assign them to the various structures
    [opts, is_valid] = assign_field(opts, field, value);
    if (~is_valid)

      [myrecording, is_valid] = assign_field(myrecording, field, value);

      if (~is_valid)
        [mysimulation, is_valid] = assign_field(mysimulation, field, value);
      end
    end

    % Otherwise we notify its inexistence
    if (~is_valid)
        warning(['Property ''' field ''' does not exist. Ignoring']);
    end
  end

  return;
end
