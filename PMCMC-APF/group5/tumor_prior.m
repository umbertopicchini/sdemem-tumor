function out = tumor_prior(theta)

% Function returning the product of priors for the unknown parameters.
% Input: theta, the current vector of parameters to be estimated. 

logbeta  = theta(1);
loggamma = theta(2);
logsigmabeta = theta(3);
logsigmaerror = theta(4);

logbetaprior = normpdf(logbeta,0.7,0.6);
gamma = exp(loggamma);
gammaprior = myinvgampdf_unnorm(gamma,5,7);
sigmabeta = exp(logsigmabeta);
sigmabetaprior = myinvgampdf_unnorm(sigmabeta,4,2);
sigmaerror = exp(logsigmaerror);
sigmaerrorprior = myinvgampdf_unnorm(sigmaerror,2,1);

% the product of the specified priors (because we assume independence on priors)
prodpriors = logbetaprior*gammaprior*sigmabetaprior*sigmaerrorprior;
% now include Jacobians for exponentiated proposals 
out = prodpriors* (gamma*sigmabeta*sigmaerror);


    
    
function y = myunifpdf(x,a,b)
    if x<a || x>b
        y=0;
    else
        y = 1/(b-a);
    end    
end

function y = mynormpdf_unnorm(x,mu,sigma) % unnormalised Gaussian pdf
    if sigma>0
       y = 1/sigma * exp(-0.5*(x-mu)^2/sigma^2);
    else
        error('the standard deviation sigma is not positive')
    end
end


function y = myinvgampdf_unnorm(x,a,b)
% the unnormalised inverse-gamma pdf with generic parameters a and b. We don't need
% the (expensive) normalising factors in Metropolis-Hastings as they
% simplify in the MH ratio

if x<=0
    y = 0;
else
    y = x^(-a-1) * exp(-b/x);
end

end


function y = mytruncgauss_unnorm(x,mu,sigma,a,b)
    % the unnormalised truncated Gaussian pdf with mean mu and SD sigma
    % truncated at [a,b]. We can remove normalization constants as they
    % simplify in the MH ratio
    if x<a || x>b
        y=0;
    else
        y = exp(-0.5*(x-mu)^2/sigma^2);  % the unnormalised density
    end
        
end


end
    


