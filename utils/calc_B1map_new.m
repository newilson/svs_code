function B1map = calc_B1map_new(image1, image2, image3, mask, alpha)

% disp ('calc_b1map3d_new begin ....');

[H,W] = size(image1);
mask1 = zeros(H,W);

im = image1;
maxval = max(im(:));
m1 = im*0.0;
m1(im > 0.1*maxval) = 1.0;
mask1 = m1;


r = image2 ./(1.0+image1);
r (image2 > image1) = -r(image2>image1);
x = 0.25 * (r + sqrt(r .* r + 8.0));
r1 = acos(x)/alpha;
r1(mask1 == 0) = 0.0;

r = image3 ./(1.0+image2);
r (image3 > image2) = -r(image3>image2);
x = 0.25 * (r + sqrt(r .* r + 8.0));
r2 = acos(x)/(2*alpha);
r2(mask1 == 0) = 0.0;

r1 (r1 < 1.0) = r2(r1 < 1.0);
r1 (r1 > 1.8) = 0;



[HM,WM] = size(mask);

if (HM > H)
    B1map = mask;
    r3 = interp2(r1);
    [H1,W1] = size(r3);
    B1map(1:H1,1:W1) = r3;
else
    B1map = r1;
end

B1map = B1map .* mask;
B1map(B1map>1.8) = 0.0;

clear r r1 r2 r3 mask1 im

% disp ('calc_b1map_new end ....');
return
end

%