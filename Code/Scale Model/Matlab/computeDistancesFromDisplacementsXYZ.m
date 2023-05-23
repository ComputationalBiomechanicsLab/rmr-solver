function distances = computeDistancesFromDisplacementsXYZ(displacementMatrixXYZ)
% Compute the distance defined by spatial displacement expressed in x, y, z
% components.
% displacementMatrixXYZ must have 3*m columns.
% USAGE distances = computeDistancesFromDisplacementsXYZ(displacementMatrixXYZ)
% Ajay Seth 2015

[nt, nc] = size(displacementMatrixXYZ);

assert(mod(nc, 3)==0, 'Input displacementMatrixXYZ must have 3x columns');
nm = nc/3;

distances = zeros(nt, nm);

for I = 1:nm,
    distances(:,I) = sqrt(sum(displacementMatrixXYZ(:,3*(I-1)+1:3*(I-1)+3).^2, 2));
end

