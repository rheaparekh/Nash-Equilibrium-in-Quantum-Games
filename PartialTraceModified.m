function Xpt = PartialTraceModified(X,varargin)
  % This is a modified version of the Inbuilt PartialTrace function in
  % QETLAB. This modification allows us to calculate the PartialTrace of a
  % symbolic matrix which is needed while generate quantum games.

  lX = length(X);
 
  % set optional argument defaults: sys=2, dim=round(sqrt(length(X))), mode=-1
  [sys,dim,mode] = opt_args({ 2, round(sqrt(lX)), -1 },varargin{:});

  num_sys = length(dim);

  % allow the user to enter a single number for dim
  if(num_sys == 1)
      dim = [dim,lX/dim];
      if abs(dim(2) - round(dim(2))) >= 2*lX*eps
          error('PartialTrace:InvalidDim','If DIM is a scalar, DIM must evenly divide length(X).');
      end
      dim(2) = round(dim(2));
      num_sys = 2;
  end

  sp = issparse(X);
  isnum = isnumeric(X);
  prod_dim = prod(dim);
  prod_dim_sys = prod(dim(sys));

  % Determine which of two computation methods to use (i.e., guess which
  % method will be faster).
  if(mode == -1)
      mode = (isnum && sp && prod_dim_sys^2 <= prod_dim);
  end

  % If the matrix is sparse and the amount we are tracing over is smaller
  % than the amount remaining, just do the naive thing and manually add up
  % the blocks.
  if(mode)
      sub_sys_vec = prod_dim*ones(1,prod_dim_sys)/prod_dim_sys;
      perm = [sys,setdiff(1:num_sys,sys)];

      X = mat2cell(PermuteSystems(X,perm,dim), sub_sys_vec, sub_sys_vec); % permute the subsystems so that we just have to do the partial trace on the first (potentially larger) subsystem
      Xpt = sparse(sub_sys_vec(1),sub_sys_vec(1));
      for j = 1:prod_dim_sys
         Xpt = Xpt + X{j,j};
      end
    
   % Otherwise, do a clever trick with mat2cell or reshaping, which is almost always faster.
   else
      sub_prod = prod_dim/prod_dim_sys;
      sub_sys_vec = prod_dim*ones(1,sub_prod)/sub_prod;
      perm = [setdiff(1:num_sys,sys),sys];

      Xpt = PermuteSystems(X,perm,dim); % permute the subsystems so that we just have to do the partial trace on the second (potentially larger) subsystem

      if(isnum) % if the input is a numeric matrix, perform the partial trace operation the fastest way we know how
          Xpt = cellfun(@(x) full(trace(x)), mat2cell(Xpt, sub_sys_vec, sub_sys_vec)); % partial trace on second subsystem
          if(sp) % if input was sparse, output should be too
              Xpt = sparse(Xpt);
          end
      else % if the input is not numeric (such as a variable in a semidefinite program), do a slower method that avoids mat2cell (mat2cell doesn't like non-numeric arrays)
          Xpt = reshape(permute(reshape(Xpt,[sub_sys_vec(1),sub_prod,sub_sys_vec(1),sub_prod]),[2,4,1,3]),[sub_prod,sub_prod,sub_sys_vec(1)^2]);
          Xpt = sumnd(Xpt(:,:,1:sub_sys_vec(1)+1:sub_sys_vec(1)^2),3);
      end
  end
end

function M=sumnd(M,dim)
   s=size(M);
   M=permute(M,[setdiff(1:ndims(M),dim),dim]);
   M=reshape(M,[],s(dim));
   M=sum(M,2);
   s(dim)=1; 
   M=reshape(M,s);
end