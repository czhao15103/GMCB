
function vech = vech(sq_mat, row)

% this function computes the vech of a square matrix

lowertriangle = tril(sq_mat);
veclt = reshape(lowertriangle, [], 1);
veclt_m0 = veclt(veclt ~= 0); % remove the zeros

if (row == true)
    vech = transpose(veclt_m0);
else 
    vech = veclt_m0;

end
