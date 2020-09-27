% Input dimensions of game
dim_a = 'Enter number of available strategies for player A: ';
dim_a = input(dim_a);
dim_b = 'Enter number of available strategies for player B: ';
dim_b = input(dim_b);
size = dim_a*dim_b;

I_a = eye(dim_a);
I_b = eye(dim_b);

% declare symbolic density matrices for both players
rho = symbolic_density_matrix(dim_a, 'x');
sigma = symbolic_density_matrix(dim_b, 'y');

% For a customised game, change the below value to required POVMs
P  = RandomPOVM(size, size,1);

% Calculate the average measurement outcome for both players
% depending on the payoff parameters and POVM selected
% Here the payoffs are randomly being selected, however you can
% customise this

payoff_A = [];
payoff_B = [];

% Comment the below code out if you are customising the payoffs
for i = 1:size
   payoff_A(i) = randi(10, 1, 1);
   payoff_B(i) = randi(10, 1, 1);
end

Ha = 0;
Hb = 0;
for i = 1:size
   Ha = Ha + payoff_A(i)*P{i};
   Hb = Hb + payoff_B(i)*P{i};
end

% Calculate the linear maps for both players

% Player 1
player_1_tensor = vpa(Tensor(rho,I_b)*Hb, 5);
phi_rho = PartialTraceModified(player_1_tensor, 1); % Applying partial trace on system 1

x = get_variables(dim_a, rho);
phi_rho(x(:)) = phi_rho;

% Player 2
player_2_tensor = vpa(Tensor(I_a, sigma)*Ha, 5);
phi_sigma = PartialTraceModified(player_2_tensor, 2); % Applying partial trace on system 2

y = get_variables(dim_b, sigma);
phi_sigma(y(:)) = phi_sigma;


function mat = symbolic_density_matrix(dim, var)
   % Generate symbolic density matrix
   mat = sym(strcat(var, '%d%d'), [dim dim]);
   for i = 1:dim
      for j = 1:dim
         if (i>j)
            s = strcat(var, int2str(j), int2str(i));
            mat(i,j) = s;
         end
      end
   end
end

function variables = get_variables(dim, mat)
   % Get a list of symbolic variables declared in matrix
   variables = sym(strcat('l%d%d'), [1 dim]);
   iter = 1;
   for i = 1:dim
      for j = 1:dim
         if (i<=j)
            variables(iter) = mat(i,j);
            iter = iter + 1;
         end
      end
   end
end  