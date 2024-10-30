function k=wvnum_omvec(h,omega,g);
% first guess
 k=(omega.*omega/g)./sqrt(tanh(omega.*omega*h/g));

%Newton Raphson

   error=ones(length(h),1)';
   while any(abs(error) > .000001)
	f=omega.*omega-g*k.*(tanh(k*h));
	fp=-g*tanh(k*h)-g*(k*h)./(cosh(k*h).^2);
	kn=k-f./fp;
	error=abs(kn-k)./k;
	k=kn;

end
	
