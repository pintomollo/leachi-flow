function [inside] = in_triangle(pts, triang)

  thresh = 2*eps;

  x = pts(:,1);
  y = pts(:,2);

  triang = triang.';

  x1 = triang(1,:);
  y1 = triang(2,:);
  x2 = triang(3,:);
  y2 = triang(4,:);
  x3 = triang(5,:);
  y3 = triang(6,:);

  s1 = bsxfun(@times, bsxfun(@minus, x, x1), bsxfun(@minus, y1, y2)) - bsxfun(@times, bsxfun(@minus, y, y1), bsxfun(@minus, x1, x2));
  s2 = bsxfun(@times, bsxfun(@minus, x, x2), bsxfun(@minus, y2, y3)) - bsxfun(@times, bsxfun(@minus, y, y2), bsxfun(@minus, x2, x3));
  s3 = bsxfun(@times, bsxfun(@minus, x, x3), bsxfun(@minus, y3, y1)) - bsxfun(@times, bsxfun(@minus, y, y3), bsxfun(@minus, x3, x1));

  s1(abs(s1) < thresh) = 0;
  s2(abs(s2) < thresh) = 0;
  s3(abs(s3) < thresh) = 0;

  s2(s2 == 0) = s1(s2 == 0);
  s2(s2 == 0) = s3(s2 == 0);

  inside = ~((xor(s1<0,s2<0) & s1) | (xor(s2<0,s3<0) & s3));

  return;
end
