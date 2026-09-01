function res = ifft2call(x)
fctr = size(x,1)*size(x,2);

X = ifftshift(ifft(fftshift(x,1),[],1),1);

res = ifftshift(ifft(fftshift(X,2),[],2),2) * sqrt(fctr);



