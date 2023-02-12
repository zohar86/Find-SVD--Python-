from cmath import sqrt
import copy
from ctypes.wintypes import FLOAT
from pickle import FALSE
from textwrap import indent, wrap
import numpy as np
from numpy import array, identity, diagonal, sin, cos, arctan, fabs
import time
from numpy import loadtxt

arrIndex = []

######## mat_identity ##########
def mat_identity(U, n):

    k = 0;
    b = n;

    for k in range(n):
        U[k][k] = 1.0
        k=k+1;
        
######## SVD ##########
def find_SVD(r, c, ATA, VT, it_max): #3,2

    n=c #n
    m=r #m

    U=np.zeros((n, n), np.float32)
    size_ATA=len(ATA)
    mat_identity(U, n)

    #U=np.identity(n, np.float32)

    flag=1
    indI = 0
    indJ = 0
    it_num = 0

   
    for it_num in range(it_max): #it_num < it_max and
        if flag:
            it_num=it_num+1

            ATA_d=copy.deepcopy(ATA)
            np.fill_diagonal(ATA_d, 0, wrap=FALSE)
        
            ATA_d[np.tril_indices(ATA_d.shape[0])] = 0
            ATA_d=np.abs(ATA_d)
            max_ind = np.argmax(ATA_d)
        
            ii=indI=int(np.floor(max_ind/size_ATA))
            jj=indJ=int(max_ind%size_ATA)
        
            max_ATA=ATA_d[ii][jj]

            if(ATA[ii][ii]==ATA[jj][jj]>0.000001):
                theta=(np.pi/4)
            else:
                g=(ATA[ii][ii]-ATA[jj][jj])
                res=(2*ATA[indI][indJ])/g
                theta=np.arctan(res)*0.5

            c=np.cos(theta)
            s=np.sin(theta)
      
            tempindI=indJ
            tempindJ=indI

            tempI=ATA[ii][ii]
            tempJ=ATA[jj][jj]
            tempMaxIJ=ATA[ii][jj]
            tempMaxJI=ATA[tempindJ][tempindI]

            g=(tempI-tempJ)

            if(np.fabs(g)<0.000001):
                g=0
            
            tempII=ATA[ii][ii] = (c * c * tempI) - (2 * s * c * tempMaxIJ) + (s * s * tempJ);
            tempJJ=ATA[jj][jj] = (s * s * tempI) + (2 * s * c * tempMaxIJ) + (c * c * tempJ);
            tempMax=ATA[ii][jj] = ATA[jj][ii] = (c * c - s * s) * (tempMaxIJ)+(s * c * (g));
            #tempMax=L[ii][jj] = L[jj][ii] = (c * c - s * s) * (tempMaxIJ)+(s * c * (m));

            p1 = indI
            p2 = 0
            q1 = indJ
            q2 = 0

            x = 0 

            itr = 0;
            arrI=np.zeros(n)

            for itr in range(n):
                arrI[x] = ATA[p1][p2];
                ATA[p1][p2] = ATA[p2][p1] = c * ATA[p1][p2] - s * ATA[q1][q2]
                x = x + 1 
                p2 = p2 + 1
                q2 = q2 + 1
                itr=itr+1

            ATA[ii][ii]=tempII
            ATA[jj][jj]=tempJJ
            ATA[ii][jj]=ATA[jj][ii]=tempMax

            p1 = indI
            p2 = 0
            q1 = indJ
            q2 = 0

            x = 0 
            itr = 0

            for i in range(n):
                ATA[q1][q2] = ATA[q2][q1] = s * arrI[x] + c * ATA[q1][q2]
                x = x + 1
                q2 = q2 + 1
                itr = itr + 1
        
            ATA[ii][ii]=tempII
            ATA[jj][jj]=tempJJ
            ATA[ii][jj]=ATA[jj][ii]=tempMax

            VindI = indI;
            VindJ = indJ;
            z = 0;
            k = 0

            for z in range(n):
                  vip = U[VindI][k];
                  viq = U[VindJ][k];
                  U[VindI][k] = c * vip - s * viq;
                  U[VindJ][k] = c * viq + s * vip;
                  k=k+1
                  z=z+1;
        
            if max_ATA<0.9:
                    flag=0
         
            #print("itr_num: ", it_num)
   
    i=0
    S= np.zeros(n, np.float32)
   
    for i in range(n):
        S[i]=ATA[i][i].real
        i=i+1
    
    VT = find_VT_Matrix(A, S, U, m, n)

    return U, S, VT, it_num
###################End################################


###########find_Matrix_VT###################### 
def find_VT_Matrix(A, S, U, m, n): 

    print("start VT: ")
    #n=row, m=col
    i=0
    j=0
    k=0

    if n>m:
        stop = m

    else:
        stop = n

    VT = np.zeros((m, m), np.float32)
    
    r=0
    for i in range(stop):
        if S[i]!=0:
            temp=1/sqrt(abs(S[i]))
            tempA=np.dot(temp, A)
            while(k<m):
                x= np.array(tempA[k])
                x=x.real
                VT[j][i]=x.dot(U[i])
                j=j+1
                k=k+1
            j=0
            k=0
        r=r+1
        print(r)
    #print("U: \n", VT)
    return VT;
#################End###############################


################ main ##########################
if __name__ == "__main__": 

    A = []
    AT = []
    ATA = [] 

    num=0
    itr_max=500000

    startTime = time.time()

    #A = loadtxt("data_81.txt", delimiter=',')
    A = loadtxt("matrixA.txt", delimiter=',')

    #matrixA.txt:
    #24, 521, 23, 23, 1
    #3, 431, 5, 32, 43
    #4, 52, 53, 52, 32
    #3, 52, 23, 5, 1
    #654, 542, 64, 14, 3

  
    #print("A: \n", A)
    A = np.array(A)
    r = len(A)
 
    AT=np.transpose(A)

    ATA= np.dot(AT, A)

    #print("ATA\n", ATA)
    c = len(AT)

    VT = np.zeros((c,c), dtype=np.float32)
    segma = np.zeros(r, dtype=np.float32)
    U = np.zeros((r,r), dtype=np.float32)

    VT, segma, U, num= find_SVD(r, c, ATA, VT, itr_max)
  
    segma=np.sqrt(np.fabs((np.sort(segma)[::-1])))  
    segmaTemp=segma[0:r]
    endTime=time.time()

    print("The total time for", num,"iterations is: \n",endTime-startTime);

    print("VT: \n", VT)
    print("segma: \n", segmaTemp)
    print("U: \n", U)

    #with open("U_res.txt", 'w') as f:
    #    f.write("------------------------------------------- The U matrix ---------------------------\n")
    #    f.write(str(U) + "\n")
    #with open("Segma_res.txt", 'w') as f:
    #    f.write("-------------------------------------------The segma matrix  ---------------------------\n")
    #    f.write(str(segmaTemp))
    #with open("VT_res.txt", 'w') as f:
    #    f.write("-------------------------------------------The V matrix  ---------------------------\n")
    #    f.write(str(VT) + "\n")

    
    ####Python library####
    startTime = time.time()
    print("\n")
    print("Pyhton library result: ")
    u, s, vt =np.linalg.svd(A)
    print("VT: \n", vt)
    print("segma: \n", s)
    print("U: \n", u)
    endTime=time.time()

    print("The total time for Python library is: \n", endTime-startTime);

 #############End########################
