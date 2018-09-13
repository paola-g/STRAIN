function [index] = jaccard(A,B)
    I = length(intersect(A,B));
    if I == length(A) + length(B)
        index = 1;
    else
    index = I/(length(A) + length(B) - I);
    end
end