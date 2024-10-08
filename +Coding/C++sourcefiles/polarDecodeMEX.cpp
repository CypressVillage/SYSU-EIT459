/* Polar Decoder
                  
Author
   - Bashar Tahir, btahir@nt.tuwien.ac.at
     (c) 2018 Institute of Telecommunications, TU Wien.
     www.nt.tuwien.ac.at 

Decoding
   - Successive Cancellation (SC) based decoder
      -> SC
      -> List-SC
	  -> CRC-List-SC (external)
*/

#include "mex.h"
#include "math.h"
#include <string.h>

const double null = 98811220033.445911;
int logn, listSize, L, n, codedBlockSize;
int test;
double *channelLLRs;
double *frozenSet;

//SC
double *decodedBits;
double *stagesLLRs;
double *stagesPSums;

//List-SC
double *decodedBitsL;
double **stagesLLRsL;
double **stagesPSumsL;

void decodeSC();
double LLRrecSC(mwSize s, mwSize i);
int calcPSumSC(mwSize s, mwSize i);
int powA(int power);
double fOper(double a,double b);

void decodeListSC();
double LLRrecListSC(mwSize s, mwSize i);
int calcPSumListSC(mwSize s, mwSize i);

int mexPrintf(const char *message, ...);
void printI(int x);
void printD(double x);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{                      
    channelLLRs = mxGetPr(prhs[0]);
    frozenSet = mxGetPr(prhs[1]);
    codedBlockSize = mxGetScalar(prhs[2]);
    listSize = mxGetScalar(prhs[3]);
    n = mxGetN(prhs[0]);
   
    if (listSize == 1)
    {
        plhs[0] = mxCreateDoubleMatrix(1, (mwSize)n, mxREAL);
        decodedBits = mxGetPr(plhs[0]);
        decodeSC();
    }else {
        plhs[0] = mxCreateDoubleMatrix((mwSize)n+1, listSize, mxREAL);
        decodedBitsL = mxGetPr(plhs[0]);
        decodeListSC();
    }
}

// Successive Cancellation
void decodeSC()  
{
    logn = log2(n);
    stagesLLRs = new double[n*logn];
    stagesPSums = new double[n*logn];

    mwSize i;
    for (i=0; i<n*logn; i++){
        stagesLLRs[i] = null;
        stagesPSums[i] = null;
    }
    for (i=0; i<codedBlockSize; i++)
    {
       if (frozenSet[i] == 1){
           decodedBits[i] = 0;
       }else{   
           if (LLRrecSC(logn, i) > 0){
               decodedBits[i] = 0;
           }else{
               decodedBits[i] = 1;
           }
       }   
    }
    delete[] stagesLLRs;
    delete[] stagesPSums;
}


// List Successive Cancellation
void decodeListSC()  
{
    logn = log2(n);
    double *LLRs;
    double *PMs;
    double *sorted;
    double **pes;
    double temp, median;
    int *remFlagged;
    int nPaths = 1;
    
    mwSize i, ii, currnPaths, k, rem, sizeStage;
    sizeStage = n*logn * sizeof(double);
    
    stagesLLRsL = new double *[listSize];
    stagesPSumsL = new double *[listSize];
    LLRs = new double [listSize];
    PMs = new double [listSize]; 
    remFlagged = new int[listSize];
    sorted = new double[listSize*2]; 
    pes = new double *[2];
    pes[0] = new double [listSize];
    pes[1] = new double [listSize];
    
    for (L=0; L<listSize; L++){
        stagesLLRsL[L] = new double [n*logn];
        stagesPSumsL[L] = new double [n*logn];
    }
    for (L=0; L<listSize; L++){
        for (i=0; i<logn*n; i++){
            stagesLLRsL[L][i] = null;
            stagesPSumsL[L][i] = null;
        }
        PMs[L] = 0;
    }
    for (i=0; i<codedBlockSize; i++){
        for (L=0; L<nPaths; L++){
            LLRs[L] = LLRrecListSC(logn, i);
        }  
        if (frozenSet[i] == 1){
            // Frozen
            for (L=0; L<nPaths; L++){
                decodedBitsL[L*(n+1)+i+1] = 0;
                if (LLRs[L] <= 0){
                    PMs[L] = PMs[L] + fabs(LLRs[L]);
                }
            }  
        }else{
            // Info
            for (L=0; L<nPaths; L++){
                if (LLRs[L] > 0){
                    pes[0][L] = PMs[L];
                    pes[1][L] = PMs[L] + fabs(LLRs[L]);
                }else{
                    pes[0][L] = PMs[L] + fabs(LLRs[L]);
                    pes[1][L] = PMs[L];
                }
            }  
            if (nPaths < listSize) {   
                currnPaths = nPaths;
                for (L=0; L<currnPaths; L++){
                    // Copy path
                    memcpy(stagesLLRsL[nPaths], stagesLLRsL[L], sizeStage);
                    memcpy(stagesPSumsL[nPaths], stagesPSumsL[L], sizeStage);
                    for (ii=0; ii<codedBlockSize; ii++){
                       decodedBitsL[nPaths*(n+1)+ii+1] = decodedBitsL[L*(n+1)+ii+1];
                    }
                    // New
                    decodedBitsL[L*(n+1)+i+1] = 0;
                    PMs[L] = pes[0][L];
                    decodedBitsL[nPaths*(n+1)+i+1] = 1;
                    PMs[nPaths] = pes[1][L];  
                    nPaths = nPaths + 1;
                }   
            }else {
                // Calculate the median
                for (L=0; L<listSize; L++){
                    sorted[L] = pes[0][L];
                    sorted[L+listSize] = pes[1][L];
                }  
                for (ii=0; ii<listSize*2; ii++){
                    for (L=1; L<listSize*2; L++){
                        if (sorted[L-1] > sorted[L]){
                            temp = sorted[L];
                            sorted[L] = sorted[L-1];
                            sorted[L-1] = temp;
                        }
                    }
                }
                median = (sorted[listSize] + sorted[listSize-1])/2;
                k = 0;
                for (L=0; L<nPaths; L++){
                    if (pes[0][L] > median && pes[1][L] > median){
                        remFlagged[k] = L;
                        k = k + 1;
                    }
                }  
                k = 0;     
                for (L=0; L<nPaths; L++){
                    if (pes[0][L] > median && pes[1][L] < median){
                        decodedBitsL[L*(n+1)+i+1] = 1;
                        PMs[L] = pes[1][L];
                    }else if (pes[1][L] > median && pes[0][L] < median){
                        decodedBitsL[L*(n+1)+i+1] = 0;
                        PMs[L] = pes[0][L];
                    }else if (pes[1][L] < median && pes[0][L] < median){
                        // Copy path
                        rem = remFlagged[k];
                        k = k + 1; 
                        memcpy(stagesLLRsL[rem], stagesLLRsL[L], sizeStage);
                        memcpy(stagesPSumsL[rem], stagesPSumsL[L], sizeStage);  
                        for (ii=0; ii<=i; ii++){
                           decodedBitsL[rem*(n+1)+ii+1] = decodedBitsL[L*(n+1)+ii+1];
                        }
                        // New
                        decodedBitsL[L*(n+1)+i+1] = 0;
                        PMs[L] = pes[0][L];
                        decodedBitsL[rem*(n+1)+i+1] = 1;
                        PMs[rem] = pes[1][L];  
                    }
                }                  
            }          
        }
    }
    for (L=0; L<nPaths; L++){
        decodedBitsL[L*(n+1)] = PMs[L];
    }
    // Free memory
    for (L=0; L<listSize; L++){
        delete[] stagesLLRsL[L];
        delete[] stagesPSumsL[L];
    }
    delete[] stagesLLRsL;
    delete[] stagesPSumsL;
    delete[] LLRs;
    delete[] PMs;
    delete[] remFlagged;
    delete[] sorted;
    delete[] pes[0];
    delete[] pes[1];
    delete[] pes;                
}

double LLRrecListSC(mwSize s, mwSize i)
{
    if (stagesLLRsL[L][logn*i + s - 1] == null){
        double a,b,out;
        
        int powS = powA(logn - s);
        if (i % (powS*2) < powS)        
        {
            if (s - 1 == 0){
                a = channelLLRs[i];
                b = channelLLRs[i + powS];
            }else{       
                a = LLRrecListSC(s - 1, i);
                b = LLRrecListSC(s - 1, i + powS);
            }
            out =  fOper(a,b);  
        }else {
            if (s - 1 == 0){
                a = channelLLRs[i - powS];
                b = channelLLRs[i];
            }else{       
                a = LLRrecListSC(s - 1, i - powS);
                b = LLRrecListSC(s - 1, i);
            }

            out = (1-2*(calcPSumListSC(s, i-powS) % 2))*a + b;
        }
        stagesLLRsL[L][logn*i + s - 1] = out;  
        return out;
    }else{
        return stagesLLRsL[L][logn*i + s - 1];
    }
}

int calcPSumListSC(mwSize s, mwSize i)
{
    int powS = powA(logn - s);
    if (s == logn) {
        return decodedBitsL[L*(n+1)+i+1];
    }else{
        if (stagesPSumsL[L][logn*i + s] == null){
            int out;
            if (i % powS < powS/2){
                out = calcPSumListSC(s + 1, i) + calcPSumListSC(s + 1, i+powS/2);
            }else{
                out = calcPSumListSC(s + 1, i);
            }
            stagesPSumsL[L][logn*i + s] = out;
            return out;
        }else{
            return stagesPSumsL[L][logn*i + s];
        }  
    }
}

double LLRrecSC(mwSize s, mwSize i)
{
    if (stagesLLRs[logn*i + s - 1] == null){
        double a,b,out;
        
        int powS = powA(logn - s);
        if (i % (powS*2) < powS)        
        {
            if (s - 1 == 0){
                a = channelLLRs[i];
                b = channelLLRs[i + powS];
            }else{       
                a = LLRrecSC(s - 1, i);
                b = LLRrecSC(s - 1, i + powS);
            }
            out =  fOper(a,b);  
        }else {
            if (s - 1 == 0){
                a = channelLLRs[i - powS];
                b = channelLLRs[i];
            }else{       
                a = LLRrecSC(s - 1, i - powS);
                b = LLRrecSC(s - 1, i);
            }

            out = (1-2*(calcPSumSC(s, i-powS) % 2))*a + b;
        }
        stagesLLRs[logn*i + s - 1] = out;  
        return out;
    }else{
        return stagesLLRs[logn*i + s - 1];
    }
}

int calcPSumSC(mwSize s, mwSize i)
{
    int powS = powA(logn - s);
    if (s == logn) {
        return decodedBits[i];
    }else{
        if (stagesPSums[logn*i + s] == null){
            int out;
            if (i % powS < powS/2){
                out = calcPSumSC(s + 1, i) + calcPSumSC(s + 1, i+powS/2);
            }else{
                out = calcPSumSC(s + 1, i);
            }
            stagesPSums[logn*i + s] = out;
            return out;
        }else{
            return stagesPSums[logn*i + s];
        }  
    }
}
int powA(int power)
{
    if (power == 0){
        return 1;
    }else{
        int out = 1;
        for (int i = 0; i < power; i++){
            out *= 2;
        }
        return out;
    }
}

double fOper(double a,double b)
{
    int aS, bS;
    if (a >= 0){
        aS = 1;
    }else{
        aS = -1;
    }
    if (b >= 0){
        bS = 1;
    }else{
        bS = -1;
    }
    if (a*aS <= b*bS){
        return bS*a;
    }else{
        return aS*b;
   }
}

void printD(double x)
{
    char buffer [10];
    sprintf (buffer, "%lf\n", x);
    mexPrintf(buffer);   
    
}
void printI(int x)
{
    char buffer [10];
    sprintf (buffer, "%d\n", x);
    mexPrintf(buffer);   
    
}