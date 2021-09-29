

function DER = compSymDeriv(R,S,ParametersID,REQ,MO)



DER = {}; % This structure will contain the DERIVATIVES
          % of reaction j with respect to metabolite i.
          

              
    
% Computation of the derivatives
fprintf('%s', 'Computing symbolic derivatives... ');

          % symbolically evaluate all parameters to compute the derivative
          % of the kinetic equations
          for ips = 1:1:numel(ParametersID)
              current = ParametersID(ips);
              eval(['syms ', char(current)]);
          end 
         % symbolically evaluate all species       
          for ips2 = 1:1:numel(MO.Species)
              current2 = MO.Species(ips2).Name;
              eval(['syms ', current2]);
          end  
          
    %% compute symbolically defivative       
for j = 1:length(R)
    for i = 1:length(S)
        kineq = REQ(j);
        kineqSym = str2sym(kineq);
        eval(['syms ', S{i}])
        DER{i,j} = diff(kineqSym,S{i});          
    end
end

save('DER','DER')

end