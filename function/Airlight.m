function A = Airlight(HazeImg, wsz)
% estimating A channel by channel separately
for k = 1 : 3
    minImg = ordfilt2(double(HazeImg(:, :, k)), 1, ones(wsz), 'symmetric');
    A(k) = max(minImg(:));
end

