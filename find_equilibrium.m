% Set the loop value
total_iterations = 1500;

% Set the starting points. These starting points are randomly generated
rho_try = RandomDensityMatrix(dim_a, 1);
sigma_try = RandomDensityMatrix(dim_b, 1);

disp("The starting random density matrices are: ")
disp(rho_try);
disp(sigma_try);

% set error tolerance
epsilon = 1.0e-6;

% Set this as true to use linear update method
% Set this as false to use matrix exponential update method
linear_update_method = false;

%set weights
% For linear updates method try weights from 0.1 to 10
% For matrix exponential method try weights from 1 to 15
weight = 5;

flag_1 = false;
epsilon_2 = 0.1;

for iteration = 1:total_iterations
    linear_map_sigma = get_linear_map_value(phi_sigma, sigma_try, dim_b);
    
    if linear_update_method == true
       update = (I + weight*linear_map_sigma)*rho_try*(I + weight*linear_map_sigma);
    else
       update = expm(weight*linear_map_sigma)*rho_try*expm(weight*linear_map_sigma);
    end   
    matrix = update/(trace(update));
    
    [V,D] = eig(linear_map_sigma);
    best_response = extract_best_response(V,D, dim_a);
    
    if get_relative_error(rho_try, best_response)<epsilon
       flag_1 = true;
    else
       flag_1 = false;
    end
    rho_try = matrix;
    
    
    linear_map_rho = get_linear_map_value(phi_rho, rho_try, dim_a);
    

    if linear_update_method == true
        update = (I + weight*linear_map_rho)*sigma_try*(I + weight*linear_map_rho);
    else    
        update = expm(weight*linear_map_rho)*sigma_try*expm(weight*linear_map_rho);
    end    
    matrix = update/trace(update);
    
    [V,D] = eig(linear_map_rho);
    best_response = extract_best_response(V,D, dim_b);
   
    if get_relative_error(sigma_try, best_response)<epsilon
        if flag_1 == true
           display_results(iteration, rho_try, sigma_try)
           break;
        end   
   end
   sigma_try = matrix;
end   


function best_response = extract_best_response(V,D, dim)
    [max_num,max_idx] = max(D(:));
    best_response = [];
    val = floor(max_idx/dim)*dim;
    for i = 1:dim
       best_response(i) = V(val+1);
       val= val+1;
    end
    best_response = vpa(best_response'*best_response,10);
end

function rel_error = get_relative_error(matrix, p)
  abs_error = abs((vpa(vpa(matrix),10)-vpa(vpa(p),10)));
  rel_error = max(abs_error(:) ./ abs(p(:)));
end

function display_results(iteration, rho_try, sigma_try)
   disp("Total number of iterations: ");
   disp(iteration);

   disp("Equilibrium value of rho: ");
   disp(rho_try);

   disp("Equilibrium value of sigma: ");
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
