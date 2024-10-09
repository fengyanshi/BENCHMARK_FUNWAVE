#define PI 3.141592653589793238462643
#define blocDim 256
#define powOfTwo 4
#define timerCount 10
#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

__global__ void guts(const long int N, const long int M,const double deltaX,const double * x,double * b)
{
	__shared__ double thisBlock[blocDim];
	__shared__ long int k,l,blockPos,start;//,fullNeeded,needed;
	double FuncTemp = 0;
	long int threadPos;

	
	if(threadIdx.x==0){
		k = blockIdx.x+2;
		l = blockIdx.y;
		start = l*(N+1);
		blockPos = k+(N+1)*l;
		//fullNeeded = k+1 + (k+1)%2;
		//needed = min(fullNeeded,blocDim);
	}
	thisBlock[threadIdx.x] = 0;
	__syncthreads();

	threadPos = threadIdx.x;
	while(threadPos<=k){
		if(threadPos==0){
			FuncTemp += x[start+threadPos];
		}else if(threadPos==1){
			if(k==2){
				FuncTemp += deltaX*x[start+threadPos]/4.0;
			}else{
				FuncTemp += deltaX*x[start+threadPos]/2.0;
			}
		}else if(threadPos<(k-1)){
			FuncTemp += deltaX*x[start+threadPos];
		}else if(threadPos==(k-1)){
			FuncTemp += 3.0/4.0*deltaX*x[start+threadPos];
		}else if(threadPos==k){
			FuncTemp += deltaX*x[start+threadPos]/4.0;
		}else{
			FuncTemp += 0;
		}

		threadPos += blockDim.x;
	}
	thisBlock[threadIdx.x] = FuncTemp;
	__syncthreads();
	for(int i=blocDim/2;i>0;i=i/2){
		if(threadIdx.x<i){
			thisBlock[threadIdx.x] += thisBlock[threadIdx.x+i];
		}
		__syncthreads();
	}
	if(threadIdx.x==0){
		b[blockPos] = thisBlock[0];
	}
	
}
