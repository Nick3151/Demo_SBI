function score = dice_coeff(x,y)
if isempty(x) || isempty(y)
    score = NaN;
else
    overlap = intersect(x,y);
    score = 2*length(overlap)/(length(x)+length(y));
end