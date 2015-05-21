function test_vignette(parameters)

  files = dir(fullfile(parameters.ImageFolder, [parameters.ImageBaseName '*']));

  all_imgs = [];
  for i=1:length(files)

    img = imread(fullfile(parameters.ImageFolder, files(i).name));

    if (parameters.flag_Color==0)
      img = double(rgb2gray(uint8(img)));
    end

    if (isempty(all_imgs))
      all_imgs = img;
    else
      all_imgs = all_imgs + img;
    end
  end

  figure;imagesc(all_imgs)
end
