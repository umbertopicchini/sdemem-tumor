function out = tumor_prior(theta)

% Function returning the product of priors for the unknown parameters.
% Input: theta, the current vector of parameters to be estimated. 

logbeta  = theta(1);
logdelta = theta(2);
logalpha = theta(3);
loggamma = theta(4);
logtau = theta(5);
logsigmabeta = theta(6);
logsigmadelta = theta(7);
logsigmaalpha = theta(8);
logsigmaerror = theta(9);


logbetaprior = normpdf(logbeta,0.7,0.6);
logdeltaprior = normpdf(logdelta,0.7,0.6);  
alpha = exp(logalpha);
alphaprior = mytruncgauss_unnorm(alpha,0.6,0.2,0.01,1);
gamma = exp(loggamma);
gammaprior = myinvgampdf_unnorm(gamma,5,7);
tau = exp(logtau);
tauprior = myinvgampdf_unnorm(tau,5,7);
sigmabeta = exp(logsigmabeta);
sigmabetaprior = myinvgampdf_unnorm(sigmabeta,4,2);
sigmadelta = exp(logsigmadelta);
sigmadeltaprior = myinvgampdf_unnorm(sigmadelta,4,2);
sigmaalpha = exp(logsigmaalpha);
sigmaalphaprior = myinvgampdf_unnorm(sigmaalpha,5,1.5);
sigmaerror = exp(logsigmaerror);
sigmaerrorprior = myinvgampdf_unnorm(sigmaerror,2,1);

prodpriors = logbetaprior*logdeltaprior*alphaprior*gammaprior*tauprior*sigmabetaprior*sigmadeltaprior*sigmaalphaprior*sigmaerrorprior;
% now include Jacobians for exponentiated proposals 
out = prodpriors* (alpha * gamma * tau * sigmabeta * sigmadelta * sigmaalpha * sigmaerror);
    
    
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
    


