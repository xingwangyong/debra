function ims(data, thresholds)
% don't want to set color map and axis every time for imagesc
% Yi Zhang 20150514
data = squeeze( data );
if ~isreal(data)
    warning('Non-real data, converting to real')
    data = abs(data);
end
figure;
if nargin == 1
%     h = imagesc(data);    
    imshowtile(data);
elseif nargin == 2
%     h = imagesc(data, thresholds);
    imshowtile(data, thresholds)
end
axis image;axis off;colormap gray;
if nargin < 2
    colorbar off;
end
end


function imshowtile(im, threshold)
    if ~isreal(im)
        im = abs(im);
    end
    im = im(:,:,:,1);
    
    % Get number of images
    n_images = size(im, 3);
    
    % For prime numbers less than 10, display in one row
    if n_images > 1 && n_images < 10 && isprime(n_images)
        tiled_im = imtile(im, 'GridSize', [1, n_images]);
    else
        tiled_im = imtile(im);
    end
    
    % Display the image
    if nargin==1
        imagesc(tiled_im);
    else
        imagesc(tiled_im, threshold);
    end
    axis image off;
    colormap gray;
    colorbar;
end