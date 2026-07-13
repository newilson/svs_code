function B1map = calc_B1map(im1, im2, mask, alpha)

H = size(mask,1);
W = size(mask,2);

% size(mask)
% size(im1)
% size(im2)

H1 = size(im1,1);
W1 = size(im1,2);

if ( H == H1 && W == W1)
    image1 = im1;
    image2 = im2;
else
    image1 = imresize(im1,[H W]);
    image2 = imresize(im2,[H W]);
end

% NW - added to match auto masking from calc_B1map_new.m
maxval = max(image1(:));
mask1 = image1*0;
mask1(image1 > 0.1*maxval) = 1;

B1map = mask .* mask1;

for w = 1:W
    for h = 1:H
%         if (mask(h,w) > 0)
        if (B1map(h,w) > 0) % NW
            r = image2(h,w) / (image1(h,w)+1e-4);
            if ( image2(h,w) > (image1(h,w)+1e-4) )
                r = -r;
            end
            x = 0.25 * ( r + sqrt( (r * r) + 8.0));
            B1map(h,w) = acos(x) / alpha;
        end
    end
end

