function x=model_MLM(a,E)

P=a(end);
a=a(1:end-1);
%E=E./(p+(1-p)*E); % Uncomment if you want to convert endmembers to albedos
x=E*a;
x=(1-P)*x./(ones(size(x))-P*x);

end