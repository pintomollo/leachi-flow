function [mystruct, is_valid] = assign_field(mystruct, field, value)

  % Find the dots indicating a subfield
  tokens = regexp(field, '\.', 'split');
  tmp_struct = mystruct;
  is_valid = true;

  % Check if the whole structure exists
  for c = 1:length(tokens)
    subfield = tokens{c};
    index = regexp(subfield, '[\(\[{].*[\)\]}]', 'match');
    if (~isempty(index))
      subfield = subfield(1:end-length(index{1}));
    end

    if (~isfield(tmp_struct, subfield))
      is_valid = false;
      break;
    elseif (c < length(tokens))
      tmp_struct = tmp_struct.(subfield);
    end
  end

  % Assign the value
  if (is_valid)
    try
      eval(['mystruct.' field ' = value;']);
    catch ME
      is_valid = false;
    end
  end

  return;
end
