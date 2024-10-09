#define PI 3.141592653589793238462643
#define blocDim 512
#define powOfTwo 4
#define timerCount 10


// Device function used when evaluating the closed integral
__device__ double closedIntt(double s, double t, double A, double B, double C)
{
	double val;
	if(t==0){
		val = 0;
	}else{
		val = t/sqrt(A)*log(2*A*s+B*t+2*sqrt(A*(A*s*s+B*s*t+C*t*t)));
	}
	return val;
}

// Device function used when evaluating the closed integral
__device__ double closedInts(double s, double t, double A, double B, double C)
{
	double val;
	if(s==0){
		val = 0;
	}else{
		val = s/sqrt(C)*log(2*C*t+B*s+2*sqrt(C*(A*s*s+B*s*t+C*t*t)));
	}
	return val;
}

// Main GPU kernel
__global__ void guts(const long int N, const long int M,
	const double * zeta,const double * zetaX,const double * zetaY,
	const double * zetaHalf,const double * zetaXHalf,const double * zetaYHalf,
	const double * phi,const double * phiX,const double * phiHalf,const double * phiXHalf,const double * phiYHalf,
	const double * x,const double * y,const double * xHalf,const double * yHalf,
	const double * weiX,const double * weiY,const double * epsilon,const double * source, const bool * singType,const double Fr,
	const double n,double * Func)
{
	__shared__ double A, B, C,sUpper,sLower,tNegUpper,tNegLower,tPosUpper,tPosLower,xH,yH,thisBlock[blocDim];
	__shared__ long int k,l,blockPos;
	double xDiff, yNegDiff, yPosDiff,xDiffSq, yNegDiffSq, yPosDiffSq, zetaDiff,radiusSqNeg,
		radiusSqPos,K1,K2,S2,I1,I2,singInt,sqrtStuff,FuncTemp = 0;
	long int i,j,threadPos,applyBlockPos;

	// First thread initialises variables
	if(threadIdx.x==0){
		k = blockIdx.x;
		l = blockIdx.y;
		blockPos = k+(N-1)*l;	// Get position of the block for mesh half ponts
		xH = xHalf[k];
		yH = yHalf[l];
		A = 1 + zetaXHalf[blockPos]*zetaXHalf[blockPos];//
		B = 2*zetaXHalf[blockPos]*zetaYHalf[blockPos];	//
		C = 1 + zetaYHalf[blockPos]*zetaYHalf[blockPos];//
		sUpper = x[N-1]-xH;								//
		sLower = x[0]-xH;								// Calculate values needed for closed integral
		tNegUpper = y[M-1]-yH;							//
		tNegLower = y[0]-yH;							//
		tPosUpper = -y[M-1]-yH;							//
		tPosLower = -y[0]-yH;							//
	}

	// After initialising variables
	__syncthreads();


	// Have each thread sum over some values of the double integrals
	threadPos = threadIdx.x;
	i = threadPos%N;
	j = threadPos/N;
	// Each loop is one collocation point
	while(threadPos<(M*N)){		

		// Calculate nessesary values
		xDiff = x[i]-xH;
		xDiffSq = xDiff*xDiff;

		yNegDiff = y[j]-yH;
		yNegDiffSq = yNegDiff*yNegDiff;

		yPosDiff = y[j]+yH;
		yPosDiffSq = yPosDiff*yPosDiff;

		zetaDiff = zeta[threadPos]-zetaHalf[blockPos];
		
		radiusSqNeg = sqrt(xDiffSq+yNegDiffSq+zetaDiff*zetaDiff);
		radiusSqPos = sqrt(xDiffSq+yPosDiffSq+zetaDiff*zetaDiff);

		// First integral kernel function
		K1 = (zetaDiff-xDiff*zetaX[threadPos]-yNegDiff*zetaY[threadPos])/(radiusSqNeg*radiusSqNeg*radiusSqNeg)+
			(zetaDiff-xDiff*zetaX[threadPos]-yPosDiff*zetaY[threadPos])/(radiusSqPos*radiusSqPos*radiusSqPos);

		// Second integral kernel function
		K2 = 1/radiusSqNeg+1/radiusSqPos;


		S2 = 1/sqrt(A*xDiffSq+B*xDiff*yNegDiff+C*yNegDiffSq)+1/sqrt(A*xDiffSq-B*xDiff*yPosDiff+C*yPosDiffSq);

		// Calculate contributions to the first and second integrals
		I1 = weiX[i]*weiY[j]*(phi[threadPos]-phiHalf[blockPos]-xDiff)*K1;

		I2 = weiX[i]*weiY[j]*(zetaX[threadPos]*K2-zetaXHalf[blockPos]*S2);

		// Accumulate integral contributions
		FuncTemp -= I1;
		FuncTemp -= I2;
		
		// Update collocation point
		threadPos += blockDim.x;
		i = threadPos%N;
		j = threadPos/N;
	}

	// Calculation of the 16 parts to the closed integral split between differnt threads
	if(threadIdx.x==blockDim.x-1){
		singInt = closedIntt(sUpper,tNegUpper,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-2){
		singInt = closedInts(sUpper,tNegUpper,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-3){
		singInt = closedIntt(sLower,tNegUpper,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-4){
		singInt = closedInts(sLower,tNegUpper,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-5){
		singInt = closedIntt(sUpper,tNegLower,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-6){
		singInt = closedInts(sUpper,tNegLower,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-7){
		singInt = closedIntt(sLower,tNegLower,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-8){
		singInt = closedInts(sLower,tNegLower,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-9){
		singInt = closedIntt(sUpper,tPosUpper,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-10){
		singInt = closedInts(sUpper,tPosUpper,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-11){
		singInt = closedIntt(sLower,tPosUpper,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-12){
		singInt = closedInts(sLower,tPosUpper,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-13){
		singInt = closedIntt(sUpper,tPosLower,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-14){
		singInt = closedInts(sUpper,tPosLower,A,B,C);
		FuncTemp -= zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-15){
		singInt = closedIntt(sLower,tPosLower,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	}else if(threadIdx.x==blockDim.x-16){
		singInt = closedInts(sLower,tPosLower,A,B,C);
		FuncTemp += zetaXHalf[blockPos]*singInt;

	// Add the phi part of the BIE
	}else if(threadIdx.x==blockDim.x-17){
		FuncTemp += 2*PI*(phiHalf[blockPos]-xHalf[k]);

	
	}
	// Add the source contribution
	if(threadIdx.x==blockDim.x-18){
		sqrtStuff = sqrt((xH-source[0])*(xH-source[0])
					+(yH-source[1])*(yH-source[1])
					+(zetaHalf[blockPos]-source[2])*(zetaHalf[blockPos]-source[2]));
			if(singType[0]){
				FuncTemp += epsilon[0]/sqrtStuff;
			}else{
				FuncTemp -= epsilon[0]*(xH-source[0])/(sqrtStuff*sqrtStuff*sqrtStuff);
			}
	}

	// Add total contribution from thread to a storage vector
	thisBlock[threadIdx.x] = FuncTemp;

	// All threads finished evaluating the BIE
	__syncthreads();
   
	// Sum up all thread contributions
	for(i=blocDim/2;i>0;i=i/2){
		if(threadIdx.x<i){
			thisBlock[threadIdx.x] += thisBlock[threadIdx.x+i];
		}
		__syncthreads();
	}

	// Store complete BIE in correct loction of output vector
	if(threadIdx.x==0){
		applyBlockPos = k+(N+1)*l+2;
		Func[applyBlockPos+M*(N+1)] = thisBlock[0];
	}

	// Split the Bernoulli's equation and radion conditions between 4 blocks

	// Fist block calculates Bernoulli's equation
	if(k==0&&l==0){

		// Have each thread compute Bernoulli's equation for a mesh half point
		threadPos = threadIdx.x;
		while(threadPos<(M*(N-1))){
			i = threadPos%(N-1);
			j = threadPos/(N-1);
			applyBlockPos = i+(N+1)*j+2;

			// Bernoulli's Equation
			Func[applyBlockPos]=((1+zetaXHalf[threadPos]*zetaXHalf[threadPos])*phiYHalf[threadPos]*phiYHalf[threadPos]
							+(1+zetaYHalf[threadPos]*zetaYHalf[threadPos])*phiXHalf[threadPos]*phiXHalf[threadPos]
							-2*zetaXHalf[threadPos]*zetaYHalf[threadPos]*phiXHalf[threadPos]*phiYHalf[threadPos])/
							(2*(1+zetaXHalf[threadPos]*zetaXHalf[threadPos]+zetaYHalf[threadPos]*zetaYHalf[threadPos]))
							+zetaHalf[threadPos]/(Fr*Fr)-0.5;
		
			threadPos += blockDim.x;
		}

	//  Second block calculates the phi radiation condition
	}else if(k==0&&l==1){
		threadPos = threadIdx.x;
		while(threadPos<M){
			i = threadPos%(N-1);
			j = threadPos/(N-1);
			applyBlockPos = (N+1)*i;

			// phi radiation condition
			Func[applyBlockPos]=x[0]*(phiX[threadPos*N]-1)+n*(phi[threadPos*N]-x[0]);
			threadPos += blockDim.x;
		}

	//  Third block calculates the phiX radiation condition
	}else if(k==0&&l==2){
		threadPos = threadIdx.x;
		while(threadPos<M){
			i = threadPos%(N-1);
			j = threadPos/(N-1);
			applyBlockPos = (N+1)*i+1;

			// phiX radiation condition
			Func[applyBlockPos]=x[0]/(x[1]-x[0])*(phiX[threadPos*N+1]-phiX[threadPos*N])+n*(phiX[threadPos*N]-1);
			threadPos += blockDim.x;
		}

	//  Fourth block calculates the zeta radiation condition
	}else if(k==0&&l==3){
		threadPos = threadIdx.x;
		while(threadPos<M){
			i = threadPos%(N-1);
			j = threadPos/(N-1);
			applyBlockPos = (N+1)*i;

			// zeta radiation condition
			Func[applyBlockPos+M*(N+1)]=x[0]*zetaX[threadPos*N]+n*zeta[threadPos*N];
			threadPos += blockDim.x;
		}

	//  Fifth block calculates the zetaX radiation condition
	}else if(k==0&&l==4){
		threadPos = threadIdx.x;
		while(threadPos<M){
			i = threadPos%(N-1);
			j = threadPos/(N-1);
			applyBlockPos = (N+1)*i+1;

			// zetaX radiation condition
			Func[applyBlockPos+M*(N+1)]=x[0]/(x[1]-x[0])*(zetaX[threadPos*N+1]-zetaX[threadPos*N])+n*zetaX[threadPos*N];
			threadPos += blockDim.x;
		}
	}

}
