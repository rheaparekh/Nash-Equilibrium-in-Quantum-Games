% Input dimensions of game
prompt = 'Enter dimension of game: ';
dim = input(prompt);

% Set the loop value
total_iterations = 1500;

% Set the starting points. These starting points are randomly generated
rho_try = RandomDensityMatrix(dim,1)
sigma_try = RandomDensityMatrix(dim,1)

flag_1 = false;

% set error tolerance
epsilon = 1.0e-6;
epsilon_2 = 0.1;

% Set this as true to use linear update method
% Set this as false to use matrix exponential update method
linear_update_method = false;

%set weights
% For linear updates method try weights from 0.1 to 10
% For matrix exponential method try weights from 1 to 15
weight = 6;


for iteration = 1:total_iterations
    linear_map_sigma = get_linear_map_value(psi_sigma, sigma_try, dim);
    
    if linear_update_method == true
       p = (I + weight*linear_map_sigma)*rho_try*(I + weight*linear_map_sigma);
    else
       p = expm(weight*linear_map_sigma)*rho_try*expm(weight*linear_map_sigma);
    end   
    p = p/(trace(p));
    
    [V,D] = eig(linear_map_sigma);
    best_response = extract_best_response(V,D, dim);
    
    if get_relative_error(rho_try, best_response)<epsilon
       flag_1 = true;
    else
       flag_1 = false;
    end
    rho_try = p;
    
    
    linear_map_rho = get_linear_map_value(psi_rho, rho_try, dim);
    

    if linear_update_method == true
        p = (I + weight*linear_map_rho)*sigma_try*(I + weight*linear_map_rho);
    else    
        p = expm(weight*linear_map_rho)*sigma_try*expm(weight*linear_map_rho);
    end    
    p = p/trace(p);
    
    [V,D] = eig(linear_map_rho);
    best_response = extract_best_response(V,D, dim);
   
    if get_relative_error(sigma_try, best_response)<epsilon
        if flag_1 == true
           display_results(iteration, rho_try, sigma_try)
           break;
        end   
   end
   sigma_try = p;
end   


function p = extract_best_response(V,D, dim)
    [max_num,max_idx] = max(D(:));
    p = [];
    val = floor(max_idx/dim)*dim;
    for i = 1:dim
       p(i) = V(val+1);
       val= val+1;
    end
    p = vpa(p'*p,10);
end

function rel_error = get_relative_error(matrix, p)
  abs_error = abs((vpa(vpa(matrix),10)-vpa(vpa(p),10)));
  rel_error = max(abs_error(:) ./ abs(p(:)))
end

function display_results(iteration, rho_try, sigma_try)
   disp(iteration);
   disp(rho_try);
   disp(sigma_try);
end

function linear_map_value = get_linear_map_value(linear_map, density_matrix, dim)
   input_vec = [];
   iter = 1;
   for i = 1:dim
      for j = 1:dim
         if (i<=j)
            input_vec(iter) = density_matrix(i,j);
            iter = iter + 1;
         end
      end
   end 
   input_value = num2cell(input_vec);
   linear_map_value = linear_map(input_value{:});
end