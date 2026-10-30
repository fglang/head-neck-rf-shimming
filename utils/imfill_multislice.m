function mask3d = imfill_multislice(mask)
%IMFILL_MULTISLICE Fill holes independently in each slice of a 3-D mask.
%   MASK3D = IMFILL_MULTISLICE(MASK) applies IMFILL(..., "holes") along
%   the third dimension. MASK and MASK3D are dimensionless binary masks.

mask3d = zeros(size(mask));
for jj = 1:size(mask3d, 3)
    mask3d(:,:,jj) = imfill(mask(:,:,jj),"holes");
end
