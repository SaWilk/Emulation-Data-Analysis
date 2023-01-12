function lgc = unfind(idx, N)
% As found on https://de.mathworks.com/matlabcentral/answers/234931-from-index-list-to-logical-index-vector
% Author: Haris K.
% Posted on: 16.02.2021
    %Go from indicies into logical (for vectors only)
    lgc = false(N,1);
    lgc(idx) = true;
end
