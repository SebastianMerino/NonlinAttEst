function [Pblock,xBlock,zBlock] = getMeanBlock(P,x,z,blockParams)
dx = x(2) - x(1);
dz = z(2) - z(1);
blockSize = blockParams.blockSize;
overlap = blockParams.overlap;

% Cropping
idx = x>blockParams.xlim(1) & x<blockParams.xlim(end);
idz = z>blockParams.zlim(1) & z<blockParams.zlim(end);
x = x(idx); z = z(idz);
P = P(idz,idx,:);

% Lateral samples
wx = round(blockSize(1)*(1-overlap)/dx);  % Between windows
nx = round(blockSize(1)/dx);                 % Window size
x0 = 1:wx:length(x)-nx;
xBlock = x(1,x0+round(nx/2));
n  = length(x0);

% Axial samples
wz = round(blockSize(2)*(1-overlap)/dz); % Between windows
nz = round(blockSize(2)/dz); % Window size
z0 = 1:wz:length(z)-nz;
zBlock = z(z0 + round(nz/2));
m  = length(z0);

Pblock = zeros(m,n,size(P,3));
for ii = 1:m
    for jj = 1:n
        blockP = P(z0(ii):z0(ii)+nz-1,x0(jj):x0(jj)+nx-1,:);
        Pblock(ii,jj,:) = mean(blockP,[1 2]);
    end
end

end