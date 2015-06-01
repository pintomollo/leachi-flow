function sizes = factorize_grid(nelems)

  center = floor(sqrt(nelems));

  sizes = [1 nelems];
  for i=2:center
    if (rem(nelems, i)==0)
      sizes(end+1,:) = [i nelems/i];
    end
  end

  return;
end
