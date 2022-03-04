function [ans1, ans2, ans3] = vectorized_for_helix(input)

for i = 1 : length(input)
    [a, b, c] = helix(input(i,1), input(i,2), input(i,3), input(i,4), input(i,5), input(i,6));
    ans (i, :) = [a, b, c]; % ans is matrix consists of raduis and pitch angles for all initial pop
end

ans1 = ans(:, 1);  %ans1 = radius
ans2 = ans(:, 2);  %ans2 = pitchangle
ans3 = ans(:, 3);  %ans3 = diameter