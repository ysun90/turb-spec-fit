function GGD = ggaussian_function(F,mu,sigma,beta)

GGD = exp( -( abs( F-mu )./(sigma.*sqrt(gamma(1/beta)/gamma(3/beta))) ).^beta ) ;
GGD=GGD/sum(GGD);

end