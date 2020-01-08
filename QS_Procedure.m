tic;
fid = fopen('QS.txt','wt');

replication = 1000;     % Number of macro replications

K = 15;                 % Number of alternatives
delta_q = 3;            % Quantile indifference zone 
d = (K-1)*delta_q/K; 
rho = 0.3;              % Induced correlation (set rho=0 if common random numbers are not used)

q = 0.9;                % Target quantile 

% Mean vector in Experiment 3
%mu = [10.5 10];

% Mean and variance vectors in Experiment 4
mu = [15 14 13 12 11 10 9 8 7 6 2 3 2 2 1];  
sigma = [0 0 0 0 0 10 9 8 7 6 3 2 2 1 1];

% These show the correct selection decisions with respect to delta_q
No_of_categories = 2;
corr_S_IN = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
             0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
selection = zeros(No_of_categories,1);

fprintf(fid, 'Number of systems: %d \n', K);
fprintf(fid, 'Quantile: %.2f \n', q);
fprintf(fid, 'Quantile IZ: %.2f \n', delta_q);
fprintf(fid, 'Number of replications: %d \n', replication);
         
sample_size = zeros(replication,1);          % Sample size information for each replication

alpha = 0.05;
alpha_k = 2*alpha/(K*(K-1)); 
g_alpha_k = 2*log(1/alpha_k); 

for l= 1:replication
    fprintf('Replication number: %d \n', l);
    
    m0 = 100;            % Initial number of observations in each batch
    b0 = 10;             % Number of batches   
    n0 = m0*b0;
    
    m = m0;
    b = b0;
    f = 2;
    
    % Variables
    S_IN = ones(1,K);     
    q_H0 = zeros(K);    % Upper triangular matrix
    q_H1 = zeros(K);    % Upper triangular matrix
        
    % Initialization
    initial_step = 1;
    X = zeros(K,b,m);   % This matrix will be updated as more observations are generated.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For Experiment 3
%       A=ones(K);
%       B=tril(A,-1);
%       C=triu(A,1);
%       v=ones(1,K);
%       Rho = rho.*B+rho.*C+diag(v);
%       Z = mvnrnd(zeros(1,K), Rho, b*m);   % bm x K
%       U = normcdf(Z,0,1);                 % bm x K
%       
%       X_all = zeros(b*m,K);
%       for i=1:K
%          X_all(:,i) =  expinv(U(:,i),mu(i));
%       end
% 
%       for i=1:b
%           X(:,i,:) = transpose(X_all((i-1)*m+1:i*m,:));           
%       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For Experiment 4
      A=ones(K);
      B=tril(A,-1);
      C=triu(A,1);
      v=ones(1,K);
      Rho = rho.*B+rho.*C+diag(v);
      Z = mvnrnd(zeros(1,K), Rho, b*m);   % bm x K
      U = normcdf(Z,0,1);                 % bm x K
      
      X_all = zeros(b*m,K);
      for i=1:5
         X_all(:,i) =  expinv(U(:,i),mu(i));
      end
      for i=6:10
          X_all(:,i) = Z(:,i)*sigma(i)+mu(i);  
      end
      for i=11:15
          X_all(:,i) = gaminv(U(:,i),mu(i),sigma(i));
      end

      for i=1:b
          X(:,i,:) = transpose(X_all((i-1)*m+1:i*m,:));           
      end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    X_merge = zeros(K,b*m);
    for i=1:K   
        for v=1:b      
            X_merge(i,(v-1)*m+1:v*m) = X(i,v,:);
        end
    end   
    X_sort = sort(X,3);    
    X_merge_sort = sort(X_merge,2);
    

    X_OBM = zeros(K,b*m-m+1,m);
    for i=1:K
        for v=1:b*m-m+1
            X_OBM(i,v,:) = X_merge(i,v:v+m-1);
        end
    end
    X_OBM_sort = sort(X_OBM,3);         % K x b*m0-m0+1 x m0
       
    sample_size(l) = K*b*m;    
    
    r = b*m;
    sampling_needed = 1;
    
    while sampling_needed 
        if initial_step == 1
            V = zeros(K,K);      % Upper triangular matrix
            T = zeros(b,m);
            T_dummy = zeros(1,m);
            t = (1:m) ./ m;
            A = zeros(K,K,b);    % Upper triangular matrix
            
            for i = 1:K
                for j = 1:K
                    if i<j                         
                        first_term = (r*m/((r-m+1)*(r-m))) * sum( ( X_OBM_sort(i,:,floor(q*m)) - X_merge_sort(i,floor(r*q)) ).^2 );
                        second_term = (r*m/((r-m+1)*(r-m))) * sum( ( X_OBM_sort(j,:,floor(q*m)) - X_merge_sort(j,floor(r*q)) ).^2 );
                        third_term = (-2)*(r*m/((r-m+1)*(r-m))) * sum( (X_OBM_sort(i,:,floor(q*m)) - X_merge_sort(i,floor(r*q))) .* (X_OBM_sort(j,:,floor(q*m)) - X_merge_sort(j,floor(r*q))) );
                        V(i,j) = first_term + second_term + third_term;                                            
                    end
                end
            end                     
        end
               
        if initial_step == 2
            X_new = zeros(K,1);         
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For Experiment 3            
%              Z = mvnrnd(zeros(1,K), Rho, 1);
%              U = normcdf(Z,0,1);
% 
%              for i=1:K
%                 if S_IN(i)== 1
%                    X_new(i) =  expinv(U(i),mu(i));
%                    sample_size(l) = sample_size(l) + 1;
%                 end
%              end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
           
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For Experiment 4
              Z = mvnrnd(zeros(1,K), Rho, 1);   % 1 x K
              U = normcdf(Z,0,1);               % 1 x K

              for i=1:K
                  if S_IN(i)== 1 
                      if i <= 5              
                         X_new(i) =  expinv(U(i),mu(i));
                      elseif i >= 6 && i <= 10
                         X_new(i) =  Z(i)*sigma(i)+mu(i); 
                      elseif i >= 11                          
                         X_new(i) = gaminv(U(i),mu(i),sigma(i));
                      end
                      sample_size(l) = sample_size(l) + 1;
                  end
              end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            X_merge = [X_merge,X_new];
            X_merge_sort = sort([X_merge_sort,X_new],2); 
            
            if sqrt(r) < m0 && r == f*n0                          
                f = 2*f;      % Since r has changed, OBM formulation will automatically calculate the new b as r-m+1.
                change_in_parameters = 1;
            elseif sqrt(r) >= m0 && floor(sqrt(r))-sqrt(r)==0
                m = sqrt(r);  % Since r has changed, OBM formulation will automatically calculate the new b as r-m+1.             
                change_in_parameters = 1;
            else
                change_in_parameters = 0;
            end

            if change_in_parameters == 1   % Update the variances            
                   
                X_OBM = zeros(K,r-m+1,m);
                for i=1:K
                    for v=1:r-m+1
                        X_OBM(i,v,:) = X_merge(i,v:v+m-1);
                    end
                end
                X_OBM_sort = sort(X_OBM,3);         % K x b*m0-m0+1 x m0               

                for i = 1:K
                    for j = 1:K
                        if i<j 
                            first_term = (r*m/((r-m+1)*(r-m))) * sum( ( X_OBM_sort(i,:,floor(q*m)) - X_merge_sort(i,floor(r*q)) ).^2 );
                            second_term = (r*m/((r-m+1)*(r-m))) * sum( ( X_OBM_sort(j,:,floor(q*m)) - X_merge_sort(j,floor(r*q)) ).^2 );
                            third_term = (-2)*(r*m/((r-m+1)*(r-m))) * sum( (X_OBM_sort(i,:,floor(q*m)) - X_merge_sort(i,floor(r*q))) .* (X_OBM_sort(j,:,floor(q*m)) - X_merge_sort(j,floor(r*q))) );
                            V(i,j) = first_term + second_term + third_term;
                        end
                    end
                end 
            end
        end
       
        % Comparison
        for i = 1:K
            for j = 1:K
                if i<j && S_IN(i)==1 && S_IN(j)==1    
                    
                    X_q_diff = X_merge_sort(i,floor(r*q)) - X_merge_sort(j,floor(r*q));   % Difference of quantile point estimates, all data
                          
                    OB_q_ij = r*(-delta_q + d) + V(i,j) * g_alpha_k/(4*d);   
                    sum_diff_X_q = r*X_q_diff;

                    if sum_diff_X_q >= OB_q_ij                                            % H0: x_q_i - x_q_j >= delta_q
                        q_H0(i,j) = 1;
                        S_IN(j) = 0;
                    elseif sum_diff_X_q <= -OB_q_ij                                       % H1: x_q_i - x_q_j <= -delta_q
                        q_H1(i,j) = 1;
                        S_IN(i) = 0;
                    end
                end %if
            end %for
        end %for

        count_S_IN = sum(S_IN);

        if count_S_IN==1 
            sampling_needed = 0;     
        else
            sampling_needed = 1;
            r = r + 1; 
            initial_step = 2; % Initialization step is completed
        end     
    end % while sampling_needed
 
    for i=1:No_of_categories
        if corr_S_IN(i,:)==S_IN 
            selection(i) = selection(i) + 1;
        end 
    end
end  % for

elapsed_time = toc;

% Output
for i=1:No_of_categories
    fprintf('Correct selection count: %d \n', selection(i));
end
fprintf('Average sample size: %6.4f \n', sum(sample_size)/replication);
fprintf('Elapsed time: %.1f \n', elapsed_time);
fprintf('Standard error of sample size: %.3f \n', sqrt(var(sample_size)/replication) );

for i=1:No_of_categories
    fprintf(fid,'Correct selection count: %d \n', selection(i));
end
fprintf(fid,'Average sample size: %6.4f \n', sum(sample_size)/replication);
fprintf(fid,'Elapsed time: %.1f \n', elapsed_time);
fprintf(fid,'Standard error of sample size: %.3f \n', sqrt(var(sample_size)/replication) );

fclose(fid);
