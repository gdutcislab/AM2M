function point = reference_point(mop)
format long     
objDim = mop.od;
% p = 3 or 2
% n = objDim+p-1，
% sample_size = C(n,p) 排列组合数
switch  mop.name
    case {'WFG1','WFG2','WFG3','WFG4','WFG5','WFG6','WFG7','WFG8','WFG9'}
        p1 = 6;
        p2 = 0;
        n1 = objDim+p1-1;
        n2 = objDim+p2-1;
        sample_size1 = nchoosek(n1,p1);
        sample_size2 = nchoosek(n2,p2);
        W1 = initweight(objDim, sample_size1);
        W2 = initweight(objDim, sample_size2);
        W2= (1 - 0.5) / objDim * ones(objDim, sample_size2) + 0.5 * W2;
        point  = [W1];        
    case {'MaOP1','MaOP2','MaOP3'}
        p1 = 2;
        p2 = 2;
        n1 = objDim+p1-1;
        n2 = objDim+p2-1;
        sample_size1 = nchoosek(n1,p1);
        sample_size2 = nchoosek(n2,p2);
        W1 = initweight(objDim, sample_size1);
        W2 = initweight(objDim, sample_size2);
        W2= (1 - 0.5) / objDim * ones(objDim, sample_size2) + 0.5 * W2;
        point  = [W1,W2];
    case{'MaOP4','MaOP5','MaOP6','MaOP7','MaOP8','MaOP9','MaOP10','MaOP11','MaOP12','MaOP13','MaOP14','MaOP15'}
        p1 = 3;
        p2 = 2;
        n1 = objDim+p1-1;
        n2 = objDim+p2-1;
        sample_size1 = nchoosek(n1,p1);
        sample_size2 = nchoosek(n2,p2);
        W1 = initweight(objDim, sample_size1);
        W2 = initweight(objDim, sample_size2);
        W2= (1 - 0.5) / objDim * ones(objDim, sample_size2) + 0.5 * W2;
        point  = [W1,W2];
end
end

% This function is written by Dr. Aimin Zhou for generating any number of weight vectors

function W = initweight(objDim, N)

    U = floor(N^(1/(objDim-1)))-2;
    M = 0;
    while M<N
        U = U+1;
        M = noweight(U, 0, objDim); 
    end

    W      = zeros(objDim, M);
    C      = 0;
    V      = zeros(objDim, 1);
    [W, C] = setweight(W, C, V, U, 0, objDim, objDim);
    W      = W / (U + 0.0);

    pos     = (W < 1.0E-5);
    W(pos)  = 1.0E-5;

end

%%
function M = noweight(unit, sum, dim)

    M = 0;

    if dim == 1
        M = 1; 
        return;
    end

    for i = 0 : 1 : (unit - sum)
        M = M + noweight(unit, sum + i, dim - 1);
    end

end

%%
function [w, c] = setweight(w, c, v, unit, sum, objdim, dim)

    if dim == objdim
        v = zeros(objdim, 1);
    end

    if dim == 1
        c       = c + 1;
        v(1)    = unit - sum;
        w(:, c)  = v;
        return;
    end

    for i = 0 : 1 : (unit - sum)
        v(dim)  = i;
        [w, c]  = setweight(w, c, v, unit, sum + i, objdim, dim - 1);
    end

end


