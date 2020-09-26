% This is a brute force method to verify the best response condition
% Here, we try to verify whether two matrices are the best responses 
% to each other 

s = 500;

for c = 1:s
    k = RandomDensityMatrix(dim_b, 1);
    if vpa(sum(dot(Hb, Tensor(rho_try, k)))) > vpa(sum(dot(Hb, Tensor(rho_try, sigma_try))))
        disp("sigma is not the best response to rho");
        disp(k);
        disp(vpa(sum(dot(Hb, Tensor(rho_try, k)))))
        break
    end
end    
     

for f = 1:s
    k = RandomDensityMatrix(dim_a, 1);
    if vpa(sum(dot(Ha, Tensor(k, sigma_try)))) > vpa(sum(dot(Ha, Tensor(rho_try, sigma_try))))
        disp("rho is not the best response to sigma");
        disp(k);
        disp(vpa(sum(dot(Ha, Tensor(k, sigma_try)))))
        break
    end
end 
