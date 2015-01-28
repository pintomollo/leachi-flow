function leachi_flow(myrecording, opts)

  if (nargin == 0)
    myrecording = [];
  elseif (nargin == 1)
    opts = get_struct('options');
  end

  if (~isstruct(myrecording))
    [myrecording, opts] = inspect_recording();

    if (~isempty(myrecording))
      [myrecording, opts] = preprocess_movie(myrecording, opts);
      save(myrecording.experiment, 'myrecording', 'opts');
    end
  end

  return;
end
