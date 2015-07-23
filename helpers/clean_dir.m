function list = clean_dir(list)

  list = list(~[list.isdir]);
  for i=length(list):-1:1
    if (list(i).name(1)=='.')
      list(i) = [];
    end
  end

  return;
end
