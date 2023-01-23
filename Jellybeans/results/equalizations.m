% read the colored image
A = imread("C:\Users\vedpi\Downloads\Caprock_color.png")
imshow(A)
% convert to 8 bit
A = im2uint8(A)
% convert to grayscale
A_gray = rgb2gray(A) 
% convert to 8 bit
A_gray = im2uint8(A_gray)
% imhist
imhist(A_gray, 100)
imshow(A_gray)
% equalization
A_eq = histeq(A_gray)
imshow(A_eq)

subplot(1,2,1); 
imshow(A_eq);
axis equal;

subplot(1,2,2);
imhist(A_eq);


% make a darker image
A_gray_dark = max(A_gray - 50, 0)
imshow(A_gray_dark)
% equalize the darker image
A_gray_dark_eq = histeq(A_gray_dark)

% see
subplot(2,3,1)
imshow(A_gray)

subplot(2,3,2)
imshow(A_eq)

subplot(2,3,3)
imhist(A_eq)

subplot(2,3,4)
imshow(A_gray_dark)

subplot(2,3,5)
imshow(A_gray_dark_eq)

subplot(2,3,6)
imhist(A_gray_dark_eq)

% read in the laplacian mask
laplace_mask = imread("C:\Users\vedpi\Downloads\Cap_gray_laplacefilter.png")
[a1, a2] = size(laplace_mask)
%laplace_mask = imresize(laplace_mask, [a1, a2]);

% add to equalized grayscale
laplace_plus_eq_gray = laplace_mask + imresize(A_eq, [a1, a2])
imshow(laplace_plus_eq_gray)

subplot(1,3,1); 
imshow(A_eq);
title('Histogram Equalized Caprock Grayscale Image');

subplot(1,3,2); 
imshow(laplace_mask);
title('Laplacian Mask');

subplot(1,3,3); 
imshow(laplace_plus_eq_gray);
title('Laplacian Mask superimposed on the Equalized Grayscale Image');

% read new cap gray
new_cap_gray = imread("C:\Users\vedpi\Downloads\Cap_gray8bit.png")
new_cap_gray_mask = imread("C:\Users\vedpi\Downloads\Cap_gray8bit_laplace.png")

cap_gray_minus_mask = new_cap_gray - new_cap_gray_mask
% add to equalized grayscale

f = figure();
f.WindowState = 'maximized';
subplot(1,3,1); 
imshow(new_cap_gray);

subplot(1,3,2); 
imshow(new_cap_gray_mask);

subplot(1,3,3); 
imshow(cap_gray_minus_mask);
