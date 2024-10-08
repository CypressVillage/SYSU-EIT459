/* LDPC Decoder
                  
Author
   - Bashar Tahir, btahir@nt.tuwien.ac.at
     (c) 2018 Institute of Telecommunications, TU Wien.
     www.nt.tuwien.ac.at 

Decoding
   - Layered Belief Propagation (LBP) based Decoder
      -> Layered Sum-Product
      -> Layered Min-Sum
	  -> Layered PWL-Min-Sum (Piecewise Linear Min-Sum)

*/

#include "mex.h"
#include "math.h"

double const INF = 10000000000;
int decodingAlgorithm;
mwSize numCN, ColCN, numVN, ColVN;

double aS, bS, cor, tempX, tempY;
double *Lj;
double *Lji;
double *cnMat;
double *vnMat;
//double *outputLji;
double *outputVarNodes;

double tanhA(double x);
double atanhA(double x);
double boxc(double x, double y);
double boxc2(double x, double y);

void decodeBP(); 
int mexPrintf(const char *message, ...);
void printI(int x);
void printD(double x);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{    
    Lj = mxGetPr(prhs[0]);
    Lji = mxGetPr(prhs[1]);
    cnMat = mxGetPr(prhs[2]);
    vnMat = mxGetPr(prhs[3]);
    numCN = mxGetScalar(prhs[4]);
    ColCN = mxGetScalar(prhs[5]);
    numVN = mxGetScalar(prhs[6]);
    ColVN = mxGetScalar(prhs[7]);
    decodingAlgorithm = mxGetScalar(prhs[8]);
    
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize)numVN, mxREAL);
    outputVarNodes = mxGetPr(plhs[0]);
    
    decodeBP();
}

// Layered Belief Propagation decoding
void decodeBP()  
{
    mwSize startInd, temp, ii, kk, k, jj, ii2, jj2, teF;
    double a, b, c;
    double *contVal = new double [ColVN];
    int *currSel = new int [ColVN];
   
    switch (decodingAlgorithm){
    // Sum-Product
        case 0:
            for (jj=0; jj<numVN; jj++){
                c = 0;
                k = -1;
                for (ii=0; ii<ColVN; ii++){
                    if (vnMat[jj+numVN*ii] != 0) {
                        k = k + 1;
                        b = 1;
                        ii2 = vnMat[jj+numVN*ii] - 1;
                        for (jj2=0; jj2<ColCN; jj2++){
                            teF = ii2+numCN*jj2;
                            if (cnMat[teF] != 0) {
                                if (cnMat[teF] - 1 != jj) {
                                     b = b * tanhA(0.5*Lji[teF]);
                                }else{
                                    currSel[k] = teF;
                                }
                            }else{
                                break;
                            }    
                        }   
                        b = 2*atanhA(b);
                        contVal[k] = b;
                        c = c + b;
                    }else{
                        break;
                    }       
                }
                outputVarNodes[jj] = Lj[jj] + c;
                
                for (ii=0; ii<=k; ii++){
                    Lji[currSel[ii]] = outputVarNodes[jj] - contVal[ii];
                }     
            }
                
            break;
            
        // Min-Sum
        case 1:
            for (jj=0; jj<numVN; jj++){
                c = 0;
                k = -1;
                for (ii=0; ii<ColVN; ii++){
                    if (vnMat[jj+numVN*ii] != 0) {
                        k = k + 1;
                        a = 1;
                        b = INF;
                        ii2 = vnMat[jj+numVN*ii] - 1;
                        for (jj2=0; jj2<ColCN; jj2++){
                            teF = ii2+numCN*jj2;
                            if (cnMat[teF] != 0) {
                                if (cnMat[teF] - 1 != jj) {
                                    b = boxc(b, Lji[teF]);
                                }else{
                                    currSel[k] = teF;
                                }
                            }else{
                                break;
                            }    
                        }   
                        contVal[k] = b;
                        c = c + b;
                    }else{
                        break;
                    }       
                }
                outputVarNodes[jj] = Lj[jj] + c;
                
                for (ii=0; ii<=k; ii++){
                    Lji[currSel[ii]] = outputVarNodes[jj] - contVal[ii];
                }     
            }
            break;
            
        // PWL-Min-Sum
        case 2:
            for (jj=0; jj<numVN; jj++){
                c = 0;
                k = -1;
                ii = 0;
                while (vnMat[jj+numVN*ii] != 0){
                    k = k + 1;
                    a = 1;
                    b = INF;
                    ii2 = vnMat[jj+numVN*ii] - 1;
                    jj2 = 0;
                    teF = ii2+numCN*jj2;
                    while (cnMat[teF] != 0) {
                        if (cnMat[teF] - 1 != jj) {
                            b = boxc2(b, Lji[teF]);
                        }else{
                            currSel[k] = teF;
                        }
                        jj2 = jj2 + 1;
                        teF = ii2+numCN*jj2;
                    }
                    contVal[k] = b;
                    c = c + b;
                    ii = ii + 1;
                }
                outputVarNodes[jj] = Lj[jj] + c;
                
                for (ii=0; ii<=k; ii++){
                    Lji[currSel[ii]] = outputVarNodes[jj] - contVal[ii];
                }     
            }
            break;
    }

    delete[] contVal;
    delete[] currSel;
    
}

double tanhA(double x){
    if (x > 10){
        return 1;
    }else if (x < -10){
        return -1;
    }else{
        return tanh(x);
    }
}
double atanhA(double x){
    if (x > 0.9999999999999999){
        return 18.7;
    }else if (x < -0.9999999999999999){
        return -18.7;
    }else{
        return atanh(x);
    }
}

double boxc(double a, double b){
    if (a >= 0){
        aS = 1;
    }else{
        aS = -1;
    }
    if (b >= 0){
        bS = 1;
    }else {
        bS = -1;
    }
    if (a*aS <= b*bS){
        return bS*a;
    }else{
        return aS*b;
   }
}

double boxc2(double a, double b){
    if (a >= 0){
        aS = 1;
    }else{
        aS = -1;
    }
    if (b >= 0){
        bS = 1;
    }else {
        bS = -1;
    }
    
    if (a*aS <= b*bS){
        return bS*a + fmax(0.6931-0.25*fabs(a+b), 0) - fmax(0.6931-0.25*fabs(a-b), 0);
    }else{
        return aS*b + fmax(0.6931-0.25*fabs(a+b), 0) - fmax(0.6931-0.25*fabs(a-b), 0);
   }
    
}

void printD(double x)
{
    char buffer [256];
    sprintf (buffer, "%lf\n", x);
    mexPrintf(buffer);   
}
void printI(int x)
{
    char buffer [256];
    sprintf (buffer, "%d\n", x);
    mexPrintf(buffer);      
}