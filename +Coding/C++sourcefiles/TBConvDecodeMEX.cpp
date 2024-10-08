/* Tail-Biting Convolutional Decoder
                  
Author
   - Bashar Tahir, btahir@nt.tuwien.ac.at
     (c) 2018 Institute of Telecommunications, TU Wien.
     www.nt.tuwien.ac.at 

Decoding
   - Maximum A Posteriori (MAP) Decoder
      -> Log-MAP
      -> MAX-Log-MAP

*/

#include "mex.h"
#include "math.h"

double const INF = 99999999999999999;
int const nStates = 64;
int const nStreams = 3;
int decodingAlgorithm;

double *currentStates;
double *nextStates;
double *vOut;
int initialState;
mwSize n, K;

double *channelLLRs;
double *outputLLRs;
double **a;
double *g;
double **b;

void decodeBCJR(); 
int mexPrintf(const char *message, ...);
void printI(int x);
void printD(double x);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{     
    channelLLRs = mxGetPr(prhs[0]);
    initialState = mxGetScalar(prhs[1]);
    currentStates = mxGetPr(prhs[2]);
    nextStates = mxGetPr(prhs[3]);
    vOut = mxGetPr(prhs[4]);
    decodingAlgorithm = mxGetScalar(prhs[5]);
    
    n = mxGetN(prhs[0]);
    K = n/nStreams;
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)K, mxREAL);
    
    outputLLRs = mxGetPr(plhs[0]);
    decodeBCJR();
}

// BCJR SISO Decoding
void decodeBCJR()  
{
    int s, k, jj, cSt, nSt, nSt2;
    double gjj, a1, a2, a3, a4, b1, b2, b3, b4, temp1, temp2, P0, P1, m1, N0, N1, m2;
    
    a = new double *[K+1];
    b = new double *[K];
    g = new double [K*nStates/2];

    for (k=0; k<K; k++){
        a[k] = new double [nStates];
        b[k] = new double [nStates];
    }
    a[K] = new double [nStates];
     
    // Initialization
    for (s=0; s<nStates; s++){
        a[0][s] = -INF;
        b[K-1][s] = -INF;
    }
    a[0][initialState] = 0;
    b[K-1][initialState] = 0;
    
    
    switch(decodingAlgorithm){
// Log-MAP
    case 0:
    // Alphas
        jj = 0;   
        for (k=0; k<K; k++){
            for (s=0; s<nStates*2; s=s+4){
                cSt = currentStates[s];
                g[jj] = vOut[cSt*6]*channelLLRs[k] + vOut[cSt*6+1]*channelLLRs[k+K] + vOut[cSt*6+2]*channelLLRs[k+2*K];

                a1 = a[k][cSt] + g[jj];
                a2 = a[k][cSt + 1] - g[jj];
                a3 = a[k][cSt] - g[jj];
                a4 = a[k][cSt + 1] + g[jj];

                a[k + 1][(int) nextStates[s]] = fmax(a1, a2) + log(1+exp(-fabs(a1-a2)));
                a[k + 1][(int) nextStates[s + 1]] = fmax(a3, a4) + log(1+exp(-fabs(a3-a4)));

                jj = jj + 1;
            }
        }
    // Betas
        for (k=K-1; k>0; k--){
            jj = 0;
            for (s=0; s<nStates*2; s=s+4){
                cSt = currentStates[s];
                nSt = nextStates[s];
                nSt2 = nextStates[s+1];
                gjj = g[k*nStates/2+jj];

                b1 = b[k][nSt] + gjj;
                b2 = b[k][nSt2] - gjj;
                b3 = b[k][nSt2] + gjj;
                b4 = b[k][nSt] - gjj;

                b[k - 1][cSt] = fmax(b1, b2) + log(1+exp(-fabs(b1-b2)));
                b[k - 1][cSt + 1] = fmax(b3, b4) + log(1+exp(-fabs(b3-b4)));

                jj = jj + 1;
            }
        }
    // Total
        jj = 0;
        for (k=0; k<K; k++){
            temp1 = -INF;
            temp2 = -INF;
            for (s=0; s<nStates*2; s=s+4){
                cSt = currentStates[s];
                nSt = nextStates[s];
                nSt2 = nextStates[s+1];

                P0 = a[k][cSt] + g[jj] + b[k][nSt];
                P1 = a[k][cSt + 1] - g[jj] + b[k][nSt];

                m1 = fmax(P0, P1) + log(1+exp(-fabs(P0-P1)));
                temp1 = fmax(temp1, m1) + log(1+exp(-fabs(m1-temp1)));

                N0 = a[k][cSt] - g[jj] + b[k][nSt2];
                N1 = a[k][cSt + 1] + g[jj] + b[k][nSt2];

                m2 = fmax(N0, N1) + log(1+exp(-fabs(N0-N1)));
                temp2 = fmax(temp2, m2) + log(1+exp(-fabs(m2-temp2)));           

                jj = jj + 1;
            }
            outputLLRs[k] = temp1-temp2;
        }
        break;
        
// MAX-Log-MAP
    case 1:     
    // Alphas
        jj = 0;   
        for (k=0; k<K; k++){
            for (s=0; s<nStates*2; s=s+4){
                cSt = currentStates[s];
                g[jj] = vOut[cSt*6]*channelLLRs[k] + vOut[cSt*6+1]*channelLLRs[k+K] + vOut[cSt*6+2]*channelLLRs[k+2*K];

                a[k + 1][(int) nextStates[s]] = fmax(a[k][cSt] + g[jj], a[k][cSt + 1] - g[jj]);
                a[k + 1][(int) nextStates[s + 1]] = fmax(a[k][cSt] - g[jj], a[k][cSt + 1] + g[jj]);

                jj = jj + 1;
            }
        }
    // Betas
        for (k=K-1; k>0; k--){
            jj = 0;
            for (s=0; s<nStates*2; s=s+4){
                cSt = currentStates[s];
                nSt = nextStates[s];
                nSt2 = nextStates[s+1];
                gjj = g[k*nStates/2+jj];

                b[k - 1][cSt] = fmax(b[k][nSt] + gjj, b[k][nSt2] - gjj);
                b[k - 1][cSt + 1] = fmax(b[k][nSt2] + gjj, b[k][nSt] - gjj);

                jj = jj + 1;
            }
        }
    // Total
        jj = 0;
        for (k=0; k<K; k++){
            temp1 = -INF;
            temp2 = -INF;
            for (s=0; s<nStates*2; s=s+4){
                cSt = currentStates[s];
                nSt = nextStates[s];
                nSt2 = nextStates[s+1];

                temp1 = fmax(temp1, fmax(a[k][cSt] + g[jj] + b[k][nSt], a[k][cSt + 1] - g[jj] + b[k][nSt]));
                temp2 = fmax(temp2, fmax(a[k][cSt] - g[jj] + b[k][nSt2], a[k][cSt + 1] + g[jj] + b[k][nSt2]));           

                jj = jj + 1;
            }
            outputLLRs[k] = temp1-temp2;
        }
        break;
    }
    
    for (k=0; k<K; k++){
        delete[] a[k];
        delete[] b[k];
    }
	delete[] a[K];
    delete[] a;
    delete[] b; 
    delete[] g;
}

void printD(double x)
{
    char buffer [100];
    sprintf (buffer, "%lf\n", x);
    mexPrintf(buffer);   
    
}
void printI(int x)
{
    char buffer [100];
    sprintf (buffer, "%d\n", x);
    mexPrintf(buffer);   
    
}