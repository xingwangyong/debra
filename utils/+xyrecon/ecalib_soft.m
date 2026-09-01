function [maps, weights] = ecalib_soft( ksp, ncalib, ksize, eigThresh_im)

eigThresh_k = 0.03; % threshold of eigenvectors in k-space

if nargin < 4
    eigThresh_im = 0.8; % threshold of eigenvectors in image space
end

% Generate ESPIRiT Maps (Takes 30 secs to 1 minute)
[sx,sy,Nc] = size(ksp);
calib = crop(ksp,[ncalib,ncalib,Nc]);

% Get maps with ESPIRiT
[kernel,S] = dat2Kernel(calib,ksize);
idx = find(S >= S(1)*eigThresh_k, 1, 'last' );
[M,W] = kernelEig(kernel(:,:,:,1:idx),[sx,sy]);
maps = M(:,:,:,end);

% Weight the eigenvectors with soft-sense eigen-values
weights = W(:,:,end) ;
weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end) > eigThresh_im);
weights = -cos(pi*weights)/2 + 1/2;
end


function res = crop(x,sx,sy,sz,st)
%  res = crop(x,sx,sy)
%  crops a 2D matrix around its center.
%
%
%  res = crop(x,sx,sy,sz,st)
%  crops a 4D matrix around its center
%
%  
%  res = crop(x,[sx,sy,sz,st])
%  same as the previous example
%
%
%
%
% (c) Michael Lustig 2007

if nargin < 2
	error('must have a target size')
end

if nargin == 2
	s = sx;
end

if nargin == 3
    s = [sx,sy];
end

if nargin == 4
    s = [sx,sy,sz];
end

if nargin == 5
    s = [sx,sy,sz,st];
end

    m = size(x);
    if length(s) < length(m)
	    s = [s, ones(1,length(m)-length(s))];
    end
	
    if sum(m==s)==length(m)
	res = x;
	return;
    end

    
    for n=1:length(s)
	    idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
    end

    % this is a dirty ugly trick
    cmd = 'res = x(idx{1}';
    for n=2:length(s)
    	cmd = sprintf('%s,idx{%d}',cmd,n);
    end
    cmd = sprintf('%s);',cmd);
    eval(cmd);
end

function [kernel, S] = dat2Kernel(data, kSize)
% kernel = dat2Kernel(data, kSize,thresh)
%
% Function to perform k-space calibration step for ESPIRiT and create
% k-space kernels. Only works for 2D multi-coil images for now.  
% 
% Inputs: 
%       data - calibration data [kx,ky,coils]
%       kSize - size of kernel (for example kSize=[6,6])
%
% Outputs: 
%       kernel - k-space kernels matrix (not cropped), which correspond to
%                the basis vectors of overlapping blocks in k-space
%       S      - (Optional parameter) The singular vectors of the
%                 calibration matrix
%
%
% See also:
%           kernelEig
%
% (c) Michael Lustig 2013



[sx,sy,nc] = size(data);
imSize = [sx,sy] ;

tmp = im2row(data,kSize); [tsx,tsy,tsz] = size(tmp);
A = reshape(tmp,tsx,tsy*tsz);

[U,S,V] = svd(A,'econ');
    
kernel = reshape(V,kSize(1),kSize(2),nc,size(V,2));
S = diag(S);S = S(:);
end

function A = im2row(k, winSize)
%res = im2row(im, winSize)
% input k (nx ny nc) : image or kspace
% input winSize (bx by) : the size of the block
% output A ((nx-bx+1)*(ny-by+1),bx*by*nc)
% Yuxin Hu, April 22,2015
[nx,ny,nc] = size(k);
bx = winSize(1);
by = winSize(2);
A = zeros((nx-bx+1)*(ny-by+1),bx*by,nc);

count = 1;
    for y = 1 : by 
        for x = 1 : bx
            A(:,count,:) = reshape(k(x:nx-bx+x,y:ny-by+y,:),(nx-bx+1)*(ny-by+1),1,nc);
            % maybe we can first reshape "k" and then use circshift to
            % accelerate it.
            count = count + 1;
        end
    end
A = reshape(A,size(A,1),size(A,2)*size(A,3));


end

function [EigenVecs, EigenVals] = kernelEig(kernel, imSize)
% [eigenVecs, eigenVals] = kernelEig(kernel, imSize)
%
% Function computes the ESPIRiT step II -- eigen-value decomposition of a 
% k-space kernel in image space. Kernels should be computed with dat2Kernel
% and then cropped to keep those corresponding to the data space. 
%
% INPUTS:
%           kernel - k-space kernels computed with dat2Kernel (4D)
%           imSize - The size of the image to compute maps for [sx,sy]
%
% OUTPUTS:
%           EigenVecs - Images representing the Eigenvectors. (sx,sy,Num coils,Num coils)
%           EigenVals - Images representing the EigenValues. (sx,sy,numcoils )
%                       The last are the largest (close to 1)
%           
% 
% See Also:
%               dat2Kernel
% 
%
% (c) Michael Lustig 2010

nc = size(kernel,3);
nv = size(kernel,4);
kSize = [size(kernel,1), size(kernel,2)];

% "rotate kernel to order by maximum variance"
k = permute(kernel,[1,2,4,3]);, k =reshape(k,prod([kSize,nv]),nc);

if size(k,1) < size(k,2)
    [u,s,v] = svd(k);
else
    
    [u,s,v] = svd(k,'econ');
end

k = k*v;
kernel = reshape(k,[kSize,nv,nc]); kernel = permute(kernel,[1,2,4,3]);


KERNEL = zeros(imSize(1), imSize(2), size(kernel,3), size(kernel,4));
for n=1:size(kernel,4)
    KERNEL(:,:,:,n) = (fft2c(zpad(conj(kernel(end:-1:1,end:-1:1,:,n))*sqrt(imSize(1)*imSize(2)), ...
        [imSize(1), imSize(2), size(kernel,3)])));
end
KERNEL = KERNEL/sqrt(prod(kSize));


EigenVecs = zeros(imSize(1), imSize(2), nc, min(nc,nv));
EigenVals = zeros(imSize(1), imSize(2), min(nc,nv));

for n=1:prod(imSize)
    [x,y] = ind2sub([imSize(1),imSize(2)],n);
    mtx = squeeze(KERNEL(x,y,:,:));

    %[C,D] = eig(mtx*mtx');
    [C,D,V] = svd(mtx,'econ');
    
    ph = repmat(exp(-i*angle(C(1,:))),[size(C,1),1]);
    C = v*(C.*ph);
    D = real(diag(D));
    EigenVals(x,y,:) = D(end:-1:1);
    EigenVecs(x,y,:,:) = C(:,end:-1:1);
end
end

function res = zpad(x,sx,sy,sz,st)
%  res = zpad(x,sx,sy)
%  Zero pads a 2D matrix around its center.
%
%
%  res = zpad(x,sx,sy,sz,st)
%  Zero pads a 4D matrix around its center
%
%  
%  res = zpad(x,[sx,sy,sz,st])
%  same as the previous example
%
%
% (c) Michael Lustig 2007

if nargin < 2
	error('must have a target size')
end

if nargin == 2
	s = sx;
end

if nargin == 3
    s = [sx,sy];
end

if nargin == 4
    s = [sx,sy,sz];
end

if nargin == 5
    s = [sx,sy,sz,st];
end

    m = size(x);
    if length(m) < length(s)
	    m = [m, ones(1,length(s)-length(m))];
    end
	
    if sum(m==s)==length(m)
	res = x;
	return;
    end

    res = zeros(s);
    
    for n=1:length(s)
	    idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
    end

    % this is a dirty ugly trick
    cmd = 'res(idx{1}';
    for n=2:length(s)
    	cmd = sprintf('%s,idx{%d}',cmd,n);
    end
    cmd = sprintf('%s)=x;',cmd);
    eval(cmd);
end
function res = fft2c(x)
fctr = size(x,1)*size(x,2);
X = fftshift(fft(ifftshift(x,1),[],1),1);
res = fftshift(fft(ifftshift(X,2),[],2),2) / sqrt(fctr);
end