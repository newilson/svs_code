function res = ifft2c(x)
fctr = size(x,1)*size(x,2);
for a=1:size(x,4)
for n=1:size(x,3)
res(:,:,n,a) = sqrt(fctr)*fftshift(ifft2(ifftshift(x(:,:,n,a))));
end
end