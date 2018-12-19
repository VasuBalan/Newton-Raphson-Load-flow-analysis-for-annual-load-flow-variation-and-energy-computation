/************************************************************************
*   Load Flow for Annual Load Pattern and Apprroximate Load pattern 	*
*								        * 
*  Date of Written : 07.01.2017					        *
*  Ref Book 1 : Computer Methods in Power system Analysis      	  	*
*  Ref Book 2 : Electrical Power System  			        *
*  https://rosettacode.org/wiki/Gaussian_elimination			* 
* 									*
*  Last Modified on : 30.10.2017					*
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265
#define mat_elem(a, y, x, n) (a + ((y) * (n) + (x)))

typedef struct {
  int nbus;
  int nline;
  int ntie;
  double Tolerance;
  double BaseMVA;
  double BaseKV;
  int MIter;
}SystemDetails;
typedef struct {
  char inputfile1[20];	// Read bus69.dat file
  char output1[20];		// ybus 
  char output2[20];		// Voltage values
  char output3[20];		// Jacobian J11, J22, J33, J44, Delp and DelQ
  char output4[20];		// Matrix inversion of J11, J22, J33 and J44
  char output5[20];		// LoadFlow Results
  char inputfile2[20];		// Read RTS Annual Load Details
  char output6[20];		// Summary of the entire.
}FileDetails;
typedef struct {
  int bno;
  int bcode_i;
  int bcode_j;
  double res;
  double rea;
  double hch;
  int sta;
}LineData;
typedef struct {
  double g;
  double b;
}Admitance;
typedef struct {
  double magni;
  double theta; // in degree
}Polar;
typedef struct {
  int bcode;
  double Ps;
  double Qs;
}BusData;
typedef struct {
  double r;
  double i;
}Rect;
typedef struct {
  int hours;
  double LoadLvl;
}AnnualLoadData;
typedef struct {
  int hours;
  double LoadLevel;
  double Ploss;
  double Qloss;
  int Vmin_busno;
  double Vmin;
}Result;
typedef struct {
  double Pload;
  double Qload;
}SystemActualLoad;


// **** Generalized Functions ****

void PrintDoubleArray(double* A, int n) {
  int i;
  for (i = 0; i < n; i++) {
    printf("\n%lf", A[i]); } 
  printf("\n"); 
}
double* initializeDoubleArray(double* A, int n) {
  int i;
  for (i = 0; i < n; i++) {
    A[i] = 0.0; }
  return A; 
}
double *createDoubleArray(double *A, int n) {
  double *array;
  array = (double *) malloc (n * sizeof(double));
  initializeDoubleArray(array, n);
  return array; 
}
BusData* initializeBusData(BusData* AB, int nb) {
  int i;
  for( i = 0; i < nb; i++) {
    AB[i].bcode = 0;
    AB[i].Ps = 0;
    AB[i].Qs = 0;
  } return AB; 
}
void printBusData(BusData* bus, int nbus) {
  int i;
  printf("\n\nBus data is as follows\nBus code\t   Generation\n\tPs(MW)\t\tQs(MVAr)\n");
  for (i = 0; i < nbus; i++) {
     printf("\n%d\t%lf\t\t%lf", bus[i].bcode, bus[i].Ps, bus[i].Qs);
  } printf("\n");
}  
double **initializeAdjMatrixN(double **A, int n) {
    int i, j;
    for(i = 0; i < n; i++) {
      for(j = 0; j < n; j++) {
	A[i][j] = 0; }
    } 
    return A; 
}
double **createAdjMatrixN(double **A, int n) {
    int i, j;
    A = (double **) malloc(n * sizeof(double *) );
    if(A == NULL) {
      printf("\nError in generalized calling function - memory allocation\n");
      exit(1);
    }
    for(i = 0; i < n; i++) {
      	A[i] = (double *) malloc ( n * sizeof(double) );
	if(A[i] == NULL) {
	    printf("\nError in generalized calling function - memory allocation\n");
	    exit(1);}
    }
    A = initializeAdjMatrixN(A, n);
    return A; 
}
void printAdjMatrixN(double **A, int n) {
  if (A == NULL) {
    printf("Matrix is NULL\n");
    } else {
      	int i, j;
	printf("\nPrinting the Matrix\n-----------------\n");
	for(i = 0; i < n; i++) {
	    for(j = 0; j < n; j++) {
		printf("%lf\t ", A[i][j]); }
	    printf("\n");
	}
	printf("-----------------\n\n"); } 
}
void PrintPolarForm(Polar* A, int n) {
  int i;
  
  for(i = 0; i < n; i++) {
    printf("\n%lf<%lf", A[i].magni, A[i].theta);
  } printf("\n");
}
Polar* InitializeDoublePolarSP(Polar *A, int nb) {
  int i,j;
  for(i = 0; i < nb; i++) {
    A[i].magni = 0.0;
    A[i].theta = 0.0;
  } return A;
}
Polar* CreatePolarStructure(Polar* A, int nb) {
  int i;
  
  A = (Polar*)malloc(nb*sizeof(Polar));
  if (A == NULL)
    printf("\nMemory allocation error - voltage polar form\n");
  
  A = InitializeDoublePolarSP(A, nb);
  
  return A; 
}
FILE *FilePointerCreation(FILE* fpo2, char *fname) {
  fpo2 = fopen(fname,"w");
  if(fpo2 == NULL) {
    printf("\nError in creating %s", fname);
    exit(1);
  }
  
  return fpo2;
}
void printDoubleArray(Polar* V, int nb) {
  int i;
  for(i = 0; i < nb; i++) {
    printf("\n%lf<%lf",V[i].magni, V[i].theta); }
  printf("\n"); 
}
void printRectArray(Rect* A , int n) {
  int i;
  
  for(i = 0; i < n; i++) {
    if(A[i].i < 0)
      printf("\n%lf%lfj", A[i].r, A[i].i);
    else
      printf("\n%lf+%lfj", A[i].r, A[i].i);
  } 
}
Rect* initializeRectDoubleArray(Rect* A, int n) {
  int i;
  
  for(i = 0; i < n; i++) {
    A[i].r = 0;
    A[i].i = 0;
  }
  return A;
}
Rect* CreateRectDoubleArray(Rect* A, int n) {
  
  A = (Rect*)malloc(n*sizeof(Rect));
  if(A == NULL)
    printf("\nError in Irect creation\n");
  
  A = initializeRectDoubleArray(A, n);
  
  return A;
}
void printRectMatrix(Rect** A, int n) {
  int i, j;
  
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      if(A[i][j].r != 0 && A[i][j].i != 0) {
	if(A[i][j].i < 0)
	  printf("\n[%d,%d] =  %lf %lfj\t",i+1, j+1, A[i][j].r, A[i][j].i);
	else
	  printf("\n[%d,%d] =  %lf + %lfj\t",i+1, j+1, A[i][j].r, A[i][j].i);
      }
    }
  } 
}
Rect** initializeRectMatrix(Rect** A, int n) {
  int i, j;
  
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      A[i][j].r = 0;
      A[i][j].i = 0;
    }
  }
  
  return A; 
}
Rect** createRectPointer(Rect** A, int n){
   int i, j;
    A = (Rect **) malloc(n * sizeof(Rect *) );
    if(A == NULL) {
      printf("\nError in Rectangular pointer - memory allocation\n");
      exit(1);
    }
    for(i = 0; i < n; i++) {
      	A[i] = (Rect *) malloc ( n * sizeof(Rect) );
	if(A[i] == NULL) {
	    printf("\nError in generalized Rectangular pointer - memory allocation\n");
	    exit(1);}
    }
    A = initializeRectMatrix(A, n);
    return A;
}
void PrintSquareIntMatrix(int** A, int n) {
  int i, j;
  for(i = 0; i < n ; i++) {
    for(j = 0; j < n ; j++) {  
      if(A[i][j] != 0)
        printf("\n[%d][%d] = %d", i+1, j+1, A[i][j]);
      //printf("%d\t", A[i][j]);  
    } // printf("\n"); 
  } 
}
int** initializeSquIntMatrix(int** A, int n) {
  int i, j;
  for(i = 0; i < n ; i++) {
    for(j = 0; j < n ; j++) {  
        A[i][j] = 0;
    }
  }
  return A;
}
int **CreateSquareIntMatrix(int** A, int n) {
  int i;
  
  A  = (int **)malloc(n * sizeof(int *));
  if(A == NULL) {
    printf("\nError in memory allocation\n\n");
    exit(1);
  }  else {
    for(i = 0; i < n; i++) {
      A[i] = (int *)malloc(n * sizeof(int));
      if(A[i] == NULL) {
        printf("\nError in memory allocation\n\n");
        exit(1);
      }
    }
  }
  
  A = initializeSquIntMatrix(A, n);
  
  return A;
}
void printIntMatrixMxN(int** A, int m, int n) {
  int i, j; 
  
  for(i = 0; i < m; i++) {  // m- lines n -buses
      for(j = 0; j < n; j++) {
        if(A[i][j] != 0)
          printf("\n[%d][%d] = %d",i+1, j+1, A[i][j]); 
      } //printf("\n");
  }  
}
int** intializeIntMatrixMxN(int** A, int m, int n) {
  int i, j; 
  
  for(i = 0; i < m; i++) {
      for(j = 0; j < n; j++) {
          A[i][j] = 0;
      }
  }
  
  printIntMatrixMxN(A, m, n);
  return A;
}
int** CreateIntMatrixMxN(int** A, int m, int n) {
  int i;
  
  A  = (int **)malloc(m * sizeof(int *));
  if(A == NULL) {
    printf("\nError in memory allocation\n\n");
    exit(1);
  }  else {
    for(i = 0; i < m; i++) {
      A[i] = (int *)malloc(n * sizeof(int));
      if(A[i] == NULL) {
        printf("\nError in memory allocation\n\n");
        exit(1);
      }
    }
  }
  
  A = intializeIntMatrixMxN(A, m, n);
  
  return A;  
}
int* initializeIntArray(int *A, int n) {
  // Initializing the array pointer
  int i;
  
  for(i = 0; i < n; i++) {
    A[i] = 0;
  }
  return A;
}
int* CreateIntArray( int *A, int n) {
  // Creating an integer array
  int i;
  
  A = (int *)malloc( n * sizeof(int) );
  if(A == NULL)
    printf("\nError in memory allocation\n\n");
 
  A = initializeIntArray(A, n);
  
  return A;
}
void printIntArray(int *A, int n) {
  int i;
  
  for(i = 0; i < n; i++) {
    printf("%d\t", A[i]);
  }
  printf("\n");
}

// **** Individual Functions ****
FileDetails ReadFileDetails(FileDetails fpd){
  // Function is to read the input and output file from offline.txt
  FILE *fp1 = fopen("Offline.txt","r");
  if(fp1 == NULL) {
    printf("\nOffline.txt doesn't exist\n");
    exit(1);
  } else {
    fscanf(fp1, "%s", fpd.inputfile1);
    fscanf(fp1, "%s", fpd.output1);
    fscanf(fp1, "%s", fpd.output2);
    fscanf(fp1, "%s", fpd.output3);
    fscanf(fp1, "%s", fpd.output4);
    fscanf(fp1, "%s", fpd.output5);
    fscanf(fp1, "%s", fpd.inputfile2);
    fscanf(fp1, "%s", fpd.output6);
  }
  fclose(fp1);
  return fpd;
}
FILE *FilePointerToRead(FILE *fp, char *fname){
  // To read the file pointer name  
  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("\nError in file openeing\n");
    exit(1);
  }
  return fp;
}
SystemDetails ReadSystemDetails(SystemDetails sys, FILE* fp1) {
  
  fscanf(fp1, "%d %d %d %lf %d %lf %lf", &sys.nbus, &sys.nline, &sys.ntie, &sys.Tolerance, &sys.MIter, &sys.BaseKV, &sys.BaseMVA);
  // printf("\n%d %d %d %lf %d %lf %lf\n", sys.nbus, sys.nline, sys.ntie, sys.Tolerance, sys.MIter, sys.BaseKV, sys.BaseMVA);
  return sys;
}
LineData *initializeLineData(LineData *AB, int nb){
  // Initializing the bus pointer
  int i;
  for( i = 0; i < nb; i++) {
	AB[i].bno = 0;  
    AB[i].bcode_i = 0;
    AB[i].bcode_j = 0;
    AB[i].res = 0;
    AB[i].rea = 0;
    AB[i].hch = 0;
	AB[i].sta =0;
  }
  return AB;
}
LineData *createLinePointer(LineData* AB, SystemDetails sys) {
  AB = (LineData*)malloc(sys.nline * sizeof(LineData));
  if(AB == NULL) {
    printf("\nError in memory allocation of bus pointer\n");
    exit(1); 
  }
  
  AB = initializeLineData(AB, sys.nline);
  return AB;
}
void printLineData(LineData *AB, int n) {
  int i;
  printf("\n\nLine data is as follows\nSNo  bcode_i  bcode_j\t\timpedance\t LCharg \tStatus ");
  for (i = 0; i < n; i++) {
    printf("\n%d\t%d\t%d\t%lf+j%lf\t%lf\t%d", AB[i].bno, AB[i].bcode_i, AB[i].bcode_j, AB[i].res, AB[i].rea, AB[i].hch, AB[i].sta);
  } printf("\n"); 
}
LineData *ReadLineData(LineData *line, FILE *fp, SystemDetails sys) {
  // Reading line data
  int i;
  
  double Zb;
  Zb = (sys.BaseKV * sys.BaseKV) / sys.BaseMVA;
    
  for (i = 0; i < sys.nline; i++) {
     fscanf(fp, "%d %d %d %lf %lf %lf %d", &line[i].bno, &line[i].bcode_i, &line[i].bcode_j, &line[i].res, &line[i].rea, &line[i].hch, &line[i].sta);
  } 
   
  for (i = 0; i < sys.nline; i++) {
    line[i].res = (line[i].res/ Zb);
    line[i].rea = (line[i].rea/ Zb);
    line[i].hch = (line[i].hch * Zb);
  } 
  
  return line;
}
LineData* copyLinePointerData(LineData* D, LineData* S, SystemDetails sys) {
  int i;
  // D - Destination  S- Source
  for(i = 0; i < sys.nline; i++) {
	 D[i].bno = S[i].bno;
	 D[i].bcode_i = S[i].bcode_i;
	 D[i].bcode_j = S[i].bcode_j;
	 D[i].res = S[i].res;
	 D[i].rea = S[i].rea;
	 D[i].hch = S[i].hch;
	 D[i].sta = S[i].sta;
  }
  return D;
}
BusData* createBusPointer(BusData* AB, SystemDetails sys) {
  AB = (BusData*)malloc(sys.nbus * sizeof(BusData));
  if(AB == NULL) {
    printf("\nError in memory allocation of bus pointer\n");
    exit(1); 
  }
  AB = initializeBusData(AB, sys.nbus);
}
BusData* ReadBusDetails(BusData* bus, FILE* fp, SystemDetails sys) {
  // Reading Bus data
  int i;
    
  for (i = 0; i < sys.nbus; i++) {
     fscanf(fp, " %d %lf %lf", &bus[i].bcode, &bus[i].Ps, &bus[i].Qs);
  }  
  return bus; 
}
BusData* copybusPointerData(BusData* D, BusData* S, SystemDetails sys) {
  // Copting bus data to tbus pointer
  int i;
  // D - Destination  S- Source  
  for (i = 0; i < sys.nbus; i++) {
     D[i].bcode = S[i].bcode;
	 D[i].Ps = S[i].Ps;
	 D[i].Qs = S[i].Qs;
  }  
  return D;	
}
Admitance **InitializeYbus(Admitance **ybus, int nb) {
  int i,j;
  for(i = 0; i < nb; i++) {
    for(j = 0; j < nb; j++) {
      ybus[i][j].g = 0.0;
      ybus[i][j].b = 0.0;
    }
  } return ybus; 
}
Admitance** CreateAdmPointer(Admitance** A, int nb) {
  int i;
   A = (Admitance**)malloc(nb*sizeof(Admitance*));
   for(i = 0; i < nb; i++) {
     A[i] = (Admitance*) malloc(nb*sizeof(Admitance)); }
   if(A == NULL) 
	printf("\nError in memory Allocation for the admitance pointer for ybus\n\n");	
   
   A = InitializeYbus(A, nb);
  return A; 
}
void printYbusMatrix(Admitance **ybus, int nb) {
  int i, j;
  printf("\nYbus matrix is as follows\n\n");
  for(i = 0; i < nb; i++) {
    for(j = 0; j < nb; j++) {
      if( (ybus[i][j].g != 0 ) && (ybus[i][j].b != 0) )
	printf("\nybus[%d][%d] = %lf+%lfj",i+1, j+1, ybus[i][j].g, ybus[i][j].b);
      if(ybus[i][j].b != 0)
	printf("\nybus[%d][%d] = %lf+%lfj",i+1, j+1, ybus[i][j].g, ybus[i][j].b);
    } 
  } printf("\n");
}
void printYbusMatrixPolar(Polar **ybs, int nb) {
  int i, j;
  printf("\nYbus matrix in polar form is as follows\n\n");
  for(i = 0; i < nb; i++) {
    for(j = 0; j < nb; j++) {
      if(ybs[i][j].magni != 0 && ybs[i][j].theta != 0)
	printf("\n%lf<%lf\t", ybs[i][j].magni, ybs[i][j].theta);
    } 
  } printf("\n"); 
}
Admitance **FormTotCharAdm(Admitance **y1, LineData* line, SystemDetails sys) {
   int i,j,k;
   double totchar;
  	
   for(i = 0; i < sys.nbus; i++) {
     for(j = 0; j < sys.nbus; j++) {
       if(i == j) { // only Diagonal Elements
	 totchar = 0.0;
	 for(k = 0; k < sys.nline; k++) {
	   if( line[k].bcode_i == i+1 )  
	     totchar =  totchar + line[k].hch; 	
	   if( line[k].bcode_j == i+1 )
	     totchar =  totchar + line[k].hch; 	
	   // printf("\nK = %d\ti=%d\ttotchar = %lf", k, i+1, totchar);
	   } y1[i][j].g = 0.0; y1[i][j].b = totchar; 
	  }	
        } }
    // printYbusMatrix(y1, nb);     
   return y1; 
}
Admitance **FormYbusElem(Admitance **y2, LineData* line, SystemDetails sys) {
  int i,j,k;
  for(i = 0; i < sys.nline; i++) {  
    if(line[i].sta != 0) {
		y2[line[i].bcode_i-1][line[i].bcode_j-1].g = line[i].res / (line[i].res*line[i].res+line[i].rea*line[i].rea);
		y2[line[i].bcode_i-1][line[i].bcode_j-1].b = - line[i].rea / (line[i].res*line[i].res+line[i].rea*line[i].rea);
    
		y2[line[i].bcode_j-1][line[i].bcode_i-1].g = y2[line[i].bcode_i-1][line[i].bcode_j-1].g;
		y2[line[i].bcode_j-1][line[i].bcode_i-1].b = y2[line[i].bcode_i-1][line[i].bcode_j-1].b;
	}
  } return y2; 
}
Admitance** Combiney1y2(Admitance ** ybus, Admitance **y1, Admitance** y2, int nb) {
  int i,j,k;
  double t1, t2;
  for(i = 0; i < nb; i++) {
     for(j = 0; j < nb; j++) {
       if(i == j) {
	 t1 = 0.0; t2 = 0.0;
	 for(k = 0; k < nb; k++) {
	   t1 = t1 + y2[i][k].g;
	   t2 = t2 + y2[i][k].b + y1[i][k].b;
	 } ybus[i][j].g = t1; ybus[i][j].b = t2;
       }
       else {
	 ybus[i][j].g = -y2[i][j].g;
	 ybus[i][j].b = -y2[i][j].b; }
     }
  } 
  return ybus; 
}
Admitance** FormationYbus(Admitance **ybus, Admitance **y1, Admitance **y2, LineData *line, SystemDetails sys) {
   int i;
   
   y1 = InitializeYbus(y1, sys.nbus);
   
   y1 = FormTotCharAdm(y1, line, sys);
   // printf("\nTotal admitance corresponding to each bus\n\n");
   // printYbusMatrix(y1, sys.nbus);	
  
   y2 = InitializeYbus(y2, sys.nbus);  
   y2 = FormYbusElem(y2, line, sys);
   // printYbusMatrix(y2, sys.nbus);	
   
   ybus = InitializeYbus(ybus, sys.nbus);
   ybus = Combiney1y2(ybus, y1, y2, sys.nbus);
   // printYbusMatrix(ybus, sys.nbus);	
   
  return ybus;
} 
Polar **InitializeYbusPolarDP(Polar **ybs, int nb) {
  int i,j;
  for(i = 0; i < nb; i++) {
    for(j = 0; j < nb; j++) {
      ybs[i][j].magni = 0.0;
      ybs[i][j].theta = 0.0;
    }
  } return ybs; 
}
Polar** ComputeMagnitudeandAngleDP(Polar **ybs, Admitance** ybus, int nb){
  int i, j;
  double val, x;
  val = 180.0 / PI;
  
  ybs = InitializeYbusPolarDP(ybs, nb);
 
  for(i = 0; i < nb; i++) {
    for(j = 0; j < nb; j++) {
      if((ybus[i][j].b == 0 && ybus[i][j].g == 0)) { // Worst case also cross checked
	ybs[i][j].magni = 0;
	ybs[i][j].theta = 0; }
      else
	ybs[i][j].magni = sqrt( (ybus[i][j].b*ybus[i][j].b) +  (ybus[i][j].g*ybus[i][j].g) );
        x = (ybus[i][j].b)/((ybus[i][j].g));
        // printf("\n%lf\t%lf\t%lf\t%lf", ybus[i][j].b, ybus[i][j].g, x, (atan(x)*180)/PI);
        if(i != j)
	  ybs[i][j].theta = (atan(x)* val) + 180;		// Now neglect
	else
	  ybs[i][j].theta = atan(x)* val;
    } }
    
    // No line exist between two nodes - NAN has to removed from angle
    for(i = 0; i < nb; i++) {
      for(j = 0; j < nb; j++) {
	if(ybs[i][j].magni == 0)
	  ybs[i][j].theta = 0; } }
/*	  
    // Negative sign in theta included // Ref Stagg Book
    for(i = 0; i < nb; i++) {
      for(j = 0; j < nb; j++) {
	if(ybs[i][j].theta != 0)
	  ybs[i][j].theta = -ybs[i][j].theta;
      } } */
 
    
  return ybs;
}
Polar** CreatePolarYbsPointer(Polar** ybs, int nb) {
  int i, j;
  
  ybs = (Polar**)malloc(nb*sizeof(Polar*));
  if(ybs == NULL) 
	printf("\nError in memory Allocation ybus - polar form\n\n");	  
  for(i = 0; i < nb; i++) {
    ybs[i] = (Polar*) malloc(nb*sizeof(Polar));
    if(ybs[i] == NULL)
      printf("\nError in memory Allocation - ybus Polar Pointer");  
  }
  
  ybs = InitializeYbusPolarDP(ybs, nb);
  return ybs;  
}
FILE* FilePointertoWrite(FILE* fpo1, Admitance** ybus, Polar** ybs, int nb) {
   int i, j;
  
    
  fprintf(fpo1,"\nYbus matrix is as follows\n\n");
  for(i = 0; i < nb; i++) {
    for(j = 0; j < nb; j++) {
      if (ybus[i][j].b != 0 && ybus[i][j].b != 0)
	if(ybus[i][j].b < 0)
	  fprintf(fpo1, "\nybus[%d][%d] = %lf%lfj", i+1, j+1, ybus[i][j].g, ybus[i][j].b);
	else
	  fprintf(fpo1, "\nybus[%d][%d] = %lf+%lfj", i+1, j+1, ybus[i][j].g, ybus[i][j].b);      
    } 
  } fprintf(fpo1,"\n"); 
  
  fprintf(fpo1, "\nYbus matrix in polar form is as follows\n\n");
  for(i = 0; i < nb; i++) {
    for(j = 0; j < nb; j++) {
      if(ybs[i][j].magni != 0 && ybs[i][j].theta != 0)
	fprintf(fpo1, "\nybs[%d][%d] = %lf<%lf", i+1, j+1, ybs[i][j].magni, ybs[i][j].theta);
    } 
  } fprintf(fpo1, "\n"); 
      
  return fpo1;
}
double* CalculatepuRealPower(double *Ppu, BusData* bus, int nb, double BaseMVA) {
  int i;
  
  Ppu = initializeDoubleArray(Ppu, nb);
  
  for(i = 0; i < nb; i++) {
    Ppu[i] = - (bus[i].Ps / (1000 * BaseMVA));
  }
  return Ppu; 
}
double* CalculatepuReacPower(double* Qpu, BusData* bus, int nb, double BaseMVA) {
  int i;
  
  Qpu = initializeDoubleArray(Qpu, nb);
    
  for(i = 0; i < nb; i++) {
    Qpu[i] = - (bus[i].Qs / (1000 * BaseMVA));
  }
  return Qpu; 
}
Polar* CovertingRecttoPolarForm(Polar* A, Rect* Volt, int nb) {		// Theta in degree
  int i;
  double x;
    
  for(i = 0; i < nb; i++) {
    A[i].magni = sqrt( (Volt[i].r * Volt[i].r) + (Volt[i].i * Volt[i].i) );
    if(Volt[i].r == 0) {
      A[i].theta = 0; A[i].theta = 0; }
    else {
      x = (Volt[i].i / Volt[i].r);	
      A[i].theta = atan(x) * 57.29577951 ; }	// α(degrees) = α(radians) × 180° / π	
  }
    
  return A;
}
Rect* InitializeVoltage(Rect* A, int n) {
  int i;
  
    
  for(i = 0; i < n; i++) {	// Flat Start voltage and zero angle 
    A[i].r = 1.0;
    A[i].i = 0.0;
  }
  return A;
}
double* CalcRealPower(double* Pcal, Polar* V, Polar** ybs, int nb) {
  int i, n, j, k;
  double temp, temp1, angle;
    
  Pcal = initializeDoubleArray(Pcal, nb);
  
  //printf("sin45 = %lf", sin(45*0.017453292));
  
  for(k = 1; k < nb; k++) {	// neglect Slack bus
    temp1 = 0;
    for(n = 0; n < nb; n++) {
      temp = 0;
      temp = ybs[k][n].magni* V[n].magni*cos((V[k].theta - V[n].theta - ybs[k][n].theta) * 0.017453292);
      temp1 = temp1 + temp;
    } Pcal[k] = V[k].magni*temp1; 
  }
  return Pcal;
}
double* CalcReacPower(double* Qcal, Polar* V, Polar** ybs, int nb) {
  int i, n, j, k;
  double temp, temp1, angle;
    
  Qcal = initializeDoubleArray(Qcal, nb);
  
  //printf("sin45 = %lf", sin(45*0.017453292));
  
  for(k = 1; k < nb; k++) {	// neglect Slack bus
    temp1 = 0;
    for(n = 0; n < nb; n++) {
      temp = 0;
      temp = ybs[k][n].magni* V[n].magni*sin((V[k].theta - V[n].theta - ybs[k][n].theta) * 0.017453292);
      temp1 = temp1 + temp;
    } Qcal[k] = V[k].magni*temp1; 
  }
  return Qcal;
}
double** CalculateJ1Matrix(double** J1, Polar* V, Polar** ybs, int nb) {
  int i, j, k, n;
  double temp, temp1;
  double VoltProd, AngleSum, angle;
  
  
  J1 = initializeAdjMatrixN(J1, nb);
  
  for(k = 1; k < nb; k++) {
    for(n = 1; n < nb; n++) {
      temp1 = 0;
      if(k == n) { 		// Diagonal Elements
	for(j = 0; j < nb; j++) {	
	  temp = 0;
	  if( j != k) {
	   temp = (ybs[k][j].magni*V[j].magni*sin((V[k].theta - V[j].theta - ybs[k][j].theta) * 0.017453292));
	   temp1 = temp1+ temp; }
	}
	J1[k][n] =  - (V[k].magni *temp1); 
      } else 		// Non Diagonal Elements
	J1[k][n] = ( (V[k].magni * ybs[k][n].magni *V[n].magni) * sin( (V[k].theta - V[n].theta - ybs[k][n].theta) * 0.017453292) );
    }
  } 
  return J1; 
}
double** CalculateJ2Matrix(double** J2, Polar* V, Polar** ybs, int nb) {
  int i, j, k, n;
  double temp, temp1;
  double VoltProd, AngleSum, angle;
  
  
  J2 = initializeAdjMatrixN(J2, nb);
  
   for(k = 1; k < nb; k++) {
    for(n = 1; n < nb; n++) {
      temp1 = 0;
      if(k == n) { 		// Diagonal Elements
	for(j = 0; j < nb; j++) {	
	  temp = 0;
	  if( j != k) {
	   temp = ( ybs[k][j].magni * (V[j].magni*cos((V[k].theta - V[j].theta - ybs[k][j].theta) * 0.017453292) ) );
	   temp1 = temp1+ temp; }
	}
	temp = 0;
	temp = ( (2 * V[k].magni * ybs[k][k].magni ) * cos((ybs[k][k].theta) * 0.017453292) ); 
	J2[k][n] =  temp + temp1; 
      } else 
	  J2[k][n] = V[k].magni * ybs[k][n].magni * cos( (V[k].theta - V[n].theta - ybs[k][n].theta) * 0.017453292);
    }
  }
  return J2; 
}
double** CalculateJ3Matrix(double** J3, Polar* V, Polar** ybs, int nb) {
  int i, j, k, n;
  double temp, temp1;
  double VoltProd, AngleSum, angle;
  
  
  J3 = initializeAdjMatrixN(J3, nb);
  
  for(k = 1; k < nb; k++) {
    for(n = 1; n < nb; n++) {
      temp1 = 0;
      if(k == n) { 		// Diagonal Elements
	for(j = 0; j < nb; j++) {	
	  temp = 0;
	  if( j != k) {
	   temp = (ybs[k][j].magni*V[j].magni*cos((V[k].theta - V[j].theta - ybs[k][j].theta) * 0.017453292));
	   temp1 = temp1+ temp; }
	}
	J3[k][n] =  V[k].magni *temp1; 
      } else 		// Non Diagonal Elements
	J3[k][n] = - ( (V[k].magni * ybs[k][n].magni *V[n].magni) * cos( (V[k].theta - V[n].theta - ybs[k][n].theta) * 0.017453292) );
    }
  } 
  return J3; 
}
double** CalculateJ4Matrix(double** J4, Polar* V, Polar** ybs, int nb) {
  int i, j, k, n;
  double temp, temp1;
  double VoltProd, AngleSum, angle;
  
  
  J4 = initializeAdjMatrixN(J4, nb);
  
   for(k = 1; k < nb; k++) {
    for(n = 1; n < nb; n++) {
      temp1 = 0;
      if(k == n) { 		// Diagonal Elements
	for(j = 0; j < nb; j++) {	
	  temp = 0;
	  if( j != k) {
	   temp = ( ybs[k][j].magni * (V[j].magni*sin((V[k].theta - V[j].theta - ybs[k][j].theta) * 0.017453292) ) );
	   temp1 = temp1+ temp; }
	}
	temp = 0;
	temp = - ( (2 * V[k].magni * ybs[k][k].magni ) * sin((ybs[k][k].theta) * 0.017453292) ); 
	J4[k][n] =  temp + temp1; 
      } else 
	  J4[k][n] = V[k].magni * ybs[k][n].magni * sin( (V[k].theta - V[n].theta - ybs[k][n].theta) * 0.017453292);
    }
  }
  return J4; 
}
double* ChangeinRealPower(double* DelP, double* Ppu, double* Pcal, int nb) {
  int i;
  
  DelP = initializeDoubleArray(DelP, nb);
  
  for(i = 1; i < nb; i++) { // Ignore slack bus
    DelP[i] = Ppu[i] - Pcal[i];
  }
  return DelP; 
}
double* ChangeinReacPower(double* DelQ, double* Qpu, double* Qcal, int nb) {
  int i;
  
  DelQ = initializeDoubleArray(DelQ, nb);
  
  for(i = 1; i < nb; i++) { // Ignore slack bus
    DelQ[i] = Qpu[i] - Qcal[i];
  }
  return DelQ; 
}
double** CopyingJtoSubMatrix(double** A1, double** A, int nb) { 
  int i, j;
  
  A1 = initializeAdjMatrixN(A1, nb-1);
  
  for(i = 1; i < nb; i++) {
    for(j = 1; j < nb; j++) {
      A1[i-1][j-1] = A[i][j];
    }
  }
  return A1; 
}
double** FinalJacobianMatrix(double** J, double** J11, double** J22, double** J33, double** J44, int nb, int Jsize) {
   int p, q, i, j, k;
   
   J = initializeAdjMatrixN(J, Jsize);
   
   j = 0;
   
   // Row Wise copying
   for(p = 0; p < Jsize ; p++) {
     for(i = 0; i < nb-1 ; i++) {
       if(i == p) {
	 for(k = 0; k < nb-1 ; k++) 	// Copying J11 
	   J[p][k] = J11[i][k];
       }
     }
     j = nb-1;
     // printf("\nj = %d", j);
     for(i = j; i < Jsize ; i++) {
       if(i == p) {
 	 for(k = 0; k < nb-1 ; k++) 	// Copying J33 
 	   J[p][k] = J33[i-j][k];
       }
     } }
   
   // Column wise copying
   for(p = nb-1; p < Jsize ; p++) {
     for(i = 0; i < nb-1 ; i++) {
       J[i][p] = J22[i][p-(nb-1)];
     }
     
     for(i = 0; i < nb-1 ; i++) {
       J[i+(nb-1)][p] = J44[i][p-(nb-1)];
     }
   }
   
   return J;
}
double* CopyingDeltoSubMatrix(double* B1, double* B, int nb) {
  int i;
  
  B1 = initializeDoubleArray(B1, nb-1);
    
  for(i = 1; i < nb; i++) {
    B1[i-1] = B[i];
  }
  return B1;
}
double* FinalKnownVector(double* x, double* DelP1, double* DelQ1,int Jsize, int nb) {
   int i, k;
   
   x = initializeDoubleArray(x, Jsize);
   
   k = 0;
   for(i = 0; i < nb-1; i++) {
     x[k] = DelP1[i];
     k++;
   }
   
   for(i = 0; i < nb-1; i++) {
     x[k] = DelQ1[i];
     k++;
   }
   
   return x; 
}
void FilePrinting(FILE* fpo1, double** A1, int Jsize, double* Aa, double* Bb,  double* A, double* B, int nb, int iter) {
  int i, j, k;
  
  fprintf(fpo1, "\n*************************\nIteration = %d\n*************************\n", iter);
  
  if (A1 == NULL) 
       printf("Matrix is NULL\n");
  else {
    fprintf(fpo1, "\nJacobian Matrix J is as folllows\n");
    fprintf(fpo1, "\n-----------------\n");
    for(i = 0; i < Jsize-1; i++) {
      for(j = 0; j < Jsize-1; j++) {
	if(A1[i][j] != 0)
	  fprintf(fpo1, "\n[%d, %d] = %lf ", i+1, j+1, A1[i][j]); }
    } fprintf(fpo1, "\n-----------------\n");    
        
    fprintf(fpo1, "\n PCalc \t QCalc \n");
    fprintf(fpo1, "\n-----------------\n");
    for(i = 1; i < nb; i++) { 	// Neglecting Slack bus
      	fprintf(fpo1, "%lf\t%lf\n ", Aa[i] , Bb[i]); }
	fprintf(fpo1, "\n-----------------\n");    
    
    fprintf(fpo1, "\nDelP \t DelQ is as follows\n");
    fprintf(fpo1, "\n-----------------\n");
    for(i = 0; i < nb-1; i++) {
      	fprintf(fpo1, "%lf\t%lf\n ", A[i] , B[i]); }
	fprintf(fpo1, "\n-----------------\n");     
   }
}
double* CopyingJtoDoubleArray(double* A1, double** A, int Jsize) { 
  int i, j, k;
    
  A1 = initializeDoubleArray(A1, (Jsize*Jsize));
  
  k = 0;
  for(i = 0; i < Jsize; i++) {
    for(j = 0; j < Jsize; j++) {
      A1[k] = A[i][j];      
        //printf("\nJele[%d] = %lf\ti= %d, j =%d", k, A[i][j], i, j);
	k = k+1;
	
    }
  }
    
  return A1;
}
void SwapRow(double *a, double *b, int r1, int r2, int n) {
  double tmp, *p1, *p2;
  int i;
  
  if (r1 == r2) return;
    for (i = 0; i < n; i++) {
      	p1 = mat_elem(a, r1, i, n);
	p2 = mat_elem(a, r2, i, n);
	tmp = *p1, *p1 = *p2, *p2 = tmp;
    }
    tmp = b[r1], b[r1] = b[r2], b[r2] = tmp; 
}
double*  GaussEliminate(double *a, double *b, double *x, int n) {
  
#define A(y, x) (*mat_elem(a, y, x, n))
  int i, j, col, row, max_row,dia;
  double max, tmp;
  
  // PrintDoubleArray(b, n);
  // PrintDoubleArray(a, (n*n));	// Test Print
  
  x = initializeDoubleArray(x, n);
  
  for (dia = 0; dia < n; dia++) {
    max_row = dia, max = A(dia, dia);
    
    for (row = dia + 1; row < n; row++)
      if ((tmp = fabs(A(row, dia))) > max)
	max_row = row, max = tmp;
      
      SwapRow(a, b, dia, max_row, n);
    
      for (row = dia + 1; row < n; row++) {
	tmp = A(row, dia) / A(dia, dia);
	for (col = dia+1; col < n; col++)
	  A(row, col) -= tmp * A(dia, col);
	  A(row, dia) = 0;
	  b[row] -= tmp * b[dia];
	}
    }
    
    for (row = n - 1; row >= 0; row--) {
      tmp = b[row];
      for (j = n - 1; j > row; j--)
	tmp -= x[j] * A(row, j);
      	x[row] = tmp / A(row, row);
      }
#undef A
return x;

}
Polar* UpdateVoltageMagnandAngle(Polar* V, double* Sol, int Jsize, int nbus) {
  int i;
  
  for(i = 1; i < nbus; i++) {
    V[i].theta = V[i].theta + (Sol[i-1]*57.29577951);		// DelDel update // Radian is changed to degree 
    V[i].magni = V[i].magni + Sol[(nbus-2)+i];	// Del V update
  }
  
  return V;
}
FILE* FilePrintingMatrixInversion(FILE *fpo4,  double* A,  int n, int iter) {
  int i;
  
  fprintf(fpo4, "\n******************************\n\tIeration = %d\n******************************\n",iter);
 
  for(i = 0; i < n; i++) {		
    fprintf(fpo4,"\n%lf", A[i]);
  } fprintf(fpo4,"\n");  
 
  return fpo4; 
}
FILE* FilePointWriteEachIter(FILE* fpo2, Polar* V,int nbus, int iter){
  int i;
  
  fprintf(fpo2, "\n******************************\n\tIeration = %d\n******************************\n",iter);
  fprintf(fpo2, "Bus Number\t Voltage\n");
  for(i = 0; i < nbus; i++) {
    fprintf(fpo2,"\n%d\t%lf<%lf", i+1, V[i].magni, V[i].theta);        
  }
  fprintf(fpo2,"\n");
  return fpo2; 
}

// **** After Convergence of load flow ****
Rect* ConvertPoltoRect(Rect* Volt, Polar* V, int nb) {
  double x, y;
  int i;
  
  for(i = 0; i < nb; i++) {
    Volt[i].r = (V[i].magni*cos( (V[i].theta) * 0.017453292 ));
    Volt[i].i = (V[i].magni*sin( (V[i].theta) * 0.017453292 ));
  }
    
  return Volt;
}
Rect* ComputeBusInjection(Rect* A, Rect* B, Admitance** C, int nb) {
  int i, j, k;
  double real, imag; 
  Rect temp;
   
  for(j = 0; j < nb; j++) {
    real = 0; imag = 0; temp.r = 0; temp.i = 0;
    for(k = 0; k < nb; k++) {
      real = 0; imag = 0;
      real = ( (C[j][k].g * B[k].r) - (C[j][k].b * B[k].i) );		// Complex multiplication of two numbers
      imag = ( (C[j][k].g * B[k].i) + (C[j][k].b * B[k].r) );
      temp.r = temp.r + real;
      temp.i = temp.i + imag;
    } A[j].r = temp.r;
    A[j].i = temp.i;
  }
  return A;
}
Rect** ComputeCurrentLineFlow(Rect** Iflow, LineData* line, Rect* A, Admitance** B, int n) {
  int i, p, q;
     
  for(i = 0; i < n; i++) {
    p = line[i].bcode_i;
    q = line[i].bcode_j;
    Iflow[p-1][q-1].r = (  ( (A[p-1].r - A[q-1].r) * B[p-1][q-1].g) - ( ((A[p-1].i - A[q-1].i)) * B[p-1][q-1].b) );
    Iflow[p-1][q-1].i = (  ( (A[p-1].r - A[q-1].r) * B[p-1][q-1].b) + ( ((A[p-1].i - A[q-1].i)) * B[p-1][q-1].g) );
    Iflow[q-1][p-1].r = - Iflow[p-1][q-1].r;
    Iflow[q-1][p-1].i = - Iflow[p-1][q-1].i;
  }
  return Iflow;
}
double conjgate(double a) {
  double temp;
  
  temp = a;
  a = -temp;
  return a;
}
Rect** ComputeLinePowerFlow(Rect** S, LineData* line, Rect** Iflow, Rect* Volt, Admitance** ybus, SystemDetails sys) {
  int i, j, p, q;
  Rect temp;
  
 // printRectMatrix(Iflow, nb); // Testing Purpose
  
  for(i = 0; i < sys.nbus; i++) {
    for(j = 0; j < sys.nbus; j++) {
      Iflow[i][j].i = conjgate(Iflow[i][j].i); // Conjugate
    }
  }
  
 // printf("\n Inside Loop\n"); printRectMatrix(Iflow, nb); // Testing Purpose
  
  for(i = 0; i < sys.nbus; i++) {
    for(j = 0; j < sys.nbus; j++) {
      if(i != j) {
	S[i][j].r = ( ( (Volt[i].r * Iflow[i][j].r) - (Volt[i].i * Iflow[i][j].i) ) * sys.BaseMVA);	// instead of BaseMVA i put 1
	S[i][j].i = ( ( (Volt[i].r * Iflow[i][j].i) + (Volt[i].i * Iflow[i][j].r) ) * sys.BaseMVA);
      }
    }
  }
  
  return S;
}
Rect* CalculateLineLoss(Rect* LineLoss, Rect** S, SystemDetails sys, LineData* line) {
  int i, p, q;
  
  for(i = 0; i < sys.nline; i++) {
    p = line[i].bcode_i; q = line[i].bcode_j;
    LineLoss[i].r = S[p-1][q-1].r + S[q-1][p-1].r;
    LineLoss[i].i = S[p-1][q-1].i + S[q-1][p-1].i;
  }
  return LineLoss;
}
Rect* CalculateSInjection(Rect* Sinj, Rect* A, Rect* B, SystemDetails sys) {
  int i, k;
  Rect temp;
  
  // S = V*I conjgate
  for(i = 0; i < sys.nbus; i++) {
    Sinj[i].r = ( ( (A[i].r * B[i].r) - (A[i].i*conjgate(B[i].i)) ) * sys.BaseMVA );
    Sinj[i].i = ( ( (A[i].r * conjgate(B[i].i)) - (A[i].i * B[i].r) ) * sys.BaseMVA );
  }
    
  return Sinj;
}
FILE* FileptrLoadFlowResults(FILE* fpo5, Rect* Volt, Rect* Sinj, BusData* bus, LineData* line, Rect** S, Rect* LineLoss, Polar* V, SystemDetails s, 
			     AnnualLoadData* ALD, int ldl, Result* Res) {
  int i, j, busno;
  
  double SumPinj, SumQinj, SumPload, SumQload, SumPLoss, SumQLoss, SPL, SQL;
  
  SumPinj = 0; SumQinj = 0; SumPload = 0; SumQload = 0; SumPLoss = 0; SumQLoss = 0;
  double minVol;
  
  
  for(i = 0; i < s.nbus; i++) {
    SumPinj = SumPinj + Sinj[i].r;
    SumQinj = SumQinj + Sinj[i].i;
    SumPload = SumPload + bus[i].Ps;
    SumQload = SumQload + bus[i].Qs;
  }
  for(i = 0; i < s.nline; i++) {
    SumPLoss = SumPLoss + LineLoss[i].r;
    SumQLoss = SumQLoss + LineLoss[i].i;
  }
  
  SPL = fabs(SumPLoss);
  SQL = fabs(SumQLoss);
  
  Res[ldl].hours = ALD[ldl].hours;
  Res[ldl].LoadLevel = ALD[ldl].LoadLvl;
  Res[ldl].Ploss = SPL*1000;
  Res[ldl].Qloss = SQL*1000;
  
  
  minVol = V[0].magni;
  for(i = 1; i < s.nbus; i++) {
    if(V[i].magni < minVol) {
      minVol = V[i].magni;
      busno = i+1;
    }
  }
  Res[ldl].Vmin_busno = busno;
  Res[ldl].Vmin = minVol;
  
  fprintf(fpo5,"%2d\t%lf\t%lf\t%lf\t%2d\t\t%lf\n", Res[ldl].hours, Res[ldl].LoadLevel, Res[ldl].Ploss, Res[ldl].Qloss, Res[ldl].Vmin_busno, Res[ldl].Vmin);
  
  
  
  return fpo5;
}
Polar* LoadFlow(Polar* V, Polar** ybs, BusData* bus, double* Ppu, double* Qpu, Rect* Volt, double* Pcal, double* Qcal, double** J1, double** J2, 
		double** J3, double** J4, double* DelP, double* DelQ, double** J11, double**J22, double** J33, double** J44, int Jsize, double**  J, 
		double *x, double* DelP1, double* DelQ1, FILE* fpo3, double* Jele, double* Sol, FILE* fpo2, FILE* fpo4, SystemDetails sys) {
  
  int iter, i;
  double Tol;
  
  iter = 1;
  Tol = 1;
        
  while((Tol > sys.Tolerance) || (iter <= sys.MIter) ) {    // Convergence loop
      if(iter == 1) {
	Volt = InitializeVoltage(Volt, sys.nbus); 	// Voltage initialization
	// printRectArray(Volt, sys.nbus);
	
	V = CovertingRecttoPolarForm(V, Volt, sys.nbus);	// Converting Rect to Polar Form theta in degree
        // PrintPolarForm(V, sys.nbus);
      }
          
      Pcal = CalcRealPower(Pcal, V, ybs, sys.nbus);
      // printf("\nPcalculated\n"); PrintDoubleArray(Pcal, sys.nbus);
          
      Qcal = CalcReacPower(Qcal, V, ybs, sys.nbus);
      // printf("\nQcalculated\n"); PrintDoubleArray(Qcal, sys.nbus);
          
      J1 = CalculateJ1Matrix(J1, V, ybs, sys.nbus);
      // printAdjMatrixN(J1, sys.nbus);
          
      J2 = CalculateJ2Matrix(J2, V, ybs, sys.nbus);
      // printAdjMatrixN(J2, sys.nbus);
          
      J3 = CalculateJ3Matrix(J3, V, ybs, sys.nbus);
      // printAdjMatrixN(J3, sys.nbus);
          
      J4 = CalculateJ4Matrix(J4, V, ybs, sys.nbus);
      // printAdjMatrixN(J4, sys.nbus);
          
      DelP = ChangeinRealPower(DelP, Ppu, Pcal, sys.nbus);
      // printf("\nDelp"); PrintDoubleArray(DelP, sys.nbus);
          
      DelQ = ChangeinReacPower(DelQ, Qpu, Qcal, sys.nbus);
      // printf("\nDelQ"); PrintDoubleArray(DelQ, sys.nbus);    
          
      J11 = CopyingJtoSubMatrix(J11, J1, sys.nbus);
      // printAdjMatrixN(J11, sys.nbus-1);
          
      J22 = CopyingJtoSubMatrix(J22, J2, sys.nbus);
      // printAdjMatrixN(J22, sys.nbus-1);
          
      J33 = CopyingJtoSubMatrix(J33, J3, sys.nbus);
      // printAdjMatrixN(J33, sys.nbus-1);
          
      J44 = CopyingJtoSubMatrix(J44, J4, sys.nbus);
      // printAdjMatrixN(J44, sys.nbus-1);
          
      J = FinalJacobianMatrix(J, J11, J22, J33, J44, sys.nbus, Jsize);
      // printAdjMatrixN(J, Jsize);
         
      DelP1 = CopyingDeltoSubMatrix(DelP1, DelP, sys.nbus);
      // printf("\nDelP1\n"); PrintDoubleArray(DelP1, sys.nbus-1);
          
      DelQ1 = CopyingDeltoSubMatrix(DelQ1, DelQ, sys.nbus);
      // printf("\nDelQ1\n"); PrintDoubleArray(DelQ1, sys.nbus-1);
          
      x = FinalKnownVector(x, DelP1, DelQ1, Jsize, sys.nbus);
      // PrintDoubleArray(x, Jsize);
          
      // Tolerance Calculation
      Tol = x[0];
      for(i = 1; i < Jsize; i++) {
      if( x[i] > Tol )
         Tol = fabs(x[i]);
      }
      
      // printf("\nTol = %lf", Tol);
          
      // Initial testing Purpose of Jacobians and DelP and DelQ
      FilePrinting(fpo3, J, Jsize, Pcal, Qcal, DelP1, DelQ1, sys.nbus, iter);
          
      // Jacobian Inverse using Gauss Elimination
      Jele = CopyingJtoDoubleArray(Jele, J, Jsize);
      // PrintDoubleArray(Jele, (Jsize*Jsize));
          
      Sol = GaussEliminate(Jele, x, Sol, Jsize);
      fpo4 = FilePrintingMatrixInversion(fpo4, x, Jsize, iter);
      // printf("\nInverse Solution\n"); PrintDoubleArray(Sol, Jsize);
          
      // FilePrintingMatrixInversion(fpo4, Sol, Jsize, iter);
      V = UpdateVoltageMagnandAngle(V, Sol, Jsize, sys.nbus);
      // printf("\nUpdate Voltage\n"); printDoubleArray(V, sys.nbus); 
          
      fpo2 = FilePointWriteEachIter(fpo2, V, sys.nbus, iter);
      
      if(Tol < sys.Tolerance)
	break;
          
      iter = iter + 1;
    }   // Convergence Loop Ends
  return V;
}
SystemActualLoad* ComputeActualLoad(SystemActualLoad* Al, SystemDetails sys, BusData* bus, int n) {
  int i;
  double pload, qload;
  
  pload = 0; qload = 0;
  
  for(i = 0; i < sys.nbus; i++) {
    pload = bus[i].Ps + pload;
    qload = bus[i].Qs + qload;
  }
  
  Al[n-1].Pload = pload;
  Al[n-1].Qload = qload;
  
  return Al;
}

// ************ Annual Load Level Data ************
AnnualLoadData* initializeLoadData(AnnualLoadData* ALD, int n) {
  int i;
  
  for(i = 0; i < n; i++) {
    ALD[i].hours = 0;
    ALD[i].LoadLvl = 0;
  }
  return ALD;
}
AnnualLoadData* CreateAnnualLoadDataPointer(AnnualLoadData* ALD, int n) {
  ALD  = (AnnualLoadData*)malloc(n*sizeof(AnnualLoadData));
  if(ALD == NULL) {
	  printf("\nError in memory allocation - Annual Load Data\n\n");
	  exit(1);
  }
  ALD = initializeLoadData(ALD, n);  
  return ALD;
}
AnnualLoadData* ReadAnnualLoadData(AnnualLoadData* ALD, FILE* fp2, int n) {
  int i;
   
  for(i = 0; i < n; i++) {
	  fscanf(fp2, "%d %lf", &ALD[i].hours, &ALD[i].LoadLvl);
  }
  return ALD;
}
void PrintAnnualLoadData(AnnualLoadData* ALD,int n) {
  int i;
  
  for(i = 0; i < n; i++) 
    printf("\n%d %lf", ALD[i].hours, ALD[i].LoadLvl);
}
BusData* UpdatePower(BusData* bus, double aa, SystemDetails sys) {
  int i;
  for(i = 0; i < sys.nbus; i++) {
      bus[i].Ps = bus[i].Ps * aa;
      bus[i].Qs = bus[i].Qs * aa;
  }
  return bus;
}

// ************ Final Result Pointer ************
SystemActualLoad* initializeActualLoadPointer(SystemActualLoad* ALoad, int n) {
  int i;
  for(i = 0; i < n; i++) {
    ALoad[i].Pload = 0;
    ALoad[i].Qload = 0;
  }
  return ALoad;
}
SystemActualLoad* CreateSystemActLoadPointer(SystemActualLoad* ALoad, int n) {
  ALoad = (SystemActualLoad*)malloc(n*sizeof(SystemActualLoad));
  if(ALoad == NULL) {
    printf("\nError in memory Allocation- Aload Pointer");
    exit(1);
  }
  ALoad = initializeActualLoadPointer(ALoad, n);
}

Result* initializeResultPointer(Result* Res, int hrs) {
  int i;
  
  for(i = 0; i < hrs; i++) {
    Res[i].hours = 0;
    Res[i].LoadLevel = 0;
    Res[i].Ploss = 0;
    Res[i].Qloss = 0;
    Res[i].Vmin_busno = 0;
    Res[i].Vmin = 0;
  }
  return Res;
}
Result* CreateResultPointer(Result* Res, int hrs) {
  Res = (Result*)malloc(hrs*sizeof(Result));
  if(Res == NULL) {
    printf("\nError in memory Allocation - Result POinter");
    exit(1);
  }
  Res = initializeResultPointer(Res, hrs);
  return Res;
}
FILE* SummaryFileWritten(FILE* fpo6, Result* Res, int hrs, SystemActualLoad* ALoad, int n, SystemDetails s) {
  int i;
  double NetPloss, NetQloss, NetEloss, NetCost;
  
  NetPloss = 0; NetQloss = 0; NetEloss = 0; NetCost = 0;
  
  fprintf(fpo6, "\nTest System : %d bus distribution system\n", s.nbus);
  fprintf(fpo6, "\nPeak load Value are : \n\tPload = %lf kW\n\tQload = %lf kVAr\n", ALoad[n-1].Pload, ALoad[n-1].Qload);
 
  
  if(hrs == 96) {
    for(i = 0; i < hrs; i++) {
      NetPloss = Res[i].Ploss + NetPloss;
      NetQloss = Res[i].Qloss + NetQloss;
    }
    fprintf(fpo6, "\nFor 96 hours load profile : \n\n\tReal Power losses : \t%lf kW\n\tReactive Power losses : \t%lf kVArn", NetPloss, NetQloss);
    fprintf(fpo6, "\n\nFor 8760 hous (yearly) load profile : \n\n\tReal Power losses : \t%lf kW\n\tReactive Power losses : \t%lf kVAr",   (NetPloss*91.25), (NetQloss*91.25));
  }

  if(hrs == 5) {
    NetPloss = (Res[0].Ploss*15) + (Res[1].Ploss*25) + (Res[2].Ploss*20) + (Res[3].Ploss*25) + (Res[4].Ploss*11);
    NetQloss = (Res[0].Qloss*15) + (Res[1].Qloss*25) + (Res[2].Qloss*20) + (Res[3].Qloss*25) + (Res[4].Qloss*11);
    fprintf(fpo6, "\nApproximate 5 Load level\n\n");    
    fprintf(fpo6, "\n\nFor 8760 hous (yearly) load profile : \n\n\tReal Power losses : \t%lf kW\n\tReactive Power losses : \t%lf kVAr",   (NetPloss*91.25), (NetQloss*91.25)); 
    
  }
  return fpo6;
}

/* ******************* Main Program Starts ***************** */

int main(int argc, char* argv[]) {
  
  int i, Jsize, hrs, ldl, one;

  hrs = 96; one = 1;
      
  FileDetails fpd = ReadFileDetails(fpd);
  
  FILE *fp1 = FilePointerToRead(fp1,fpd.inputfile1);
  
  SystemDetails sys = ReadSystemDetails(sys, fp1);
  
  LineData *line = createLinePointer(line, sys);
  // printLineData(line, sys.nline);
  
  line = ReadLineData(line, fp1, sys);
  // printLineData(line, sys.nline);
  
  LineData *tline = createLinePointer(tline, sys);
  // printLineData(tline, sys.nline);
  
  BusData* bus = createBusPointer(bus, sys);
  // printBusData(bus, sys.nbus);
   
  bus = ReadBusDetails(bus, fp1, sys);  
  // printBusData(bus, sys.nbus);
  
  // creating temporary bus pointer
  BusData* tbus = createBusPointer(tbus, sys);	
  // printBusData(tbus, sys.nbus);
  
  //Before computing per unit power then copied to another pointer
  tbus = copybusPointerData(tbus, bus, sys);    
  // printBusData(tbus, sys.nbus);
  
  Admitance** y1 = CreateAdmPointer(y1, sys.nbus);  // Self Element Pointer
  
  Admitance** y2 = CreateAdmPointer(y2, sys.nbus);  // Mutual Elemnet Pointer
  
  Admitance** ybus = CreateAdmPointer(ybus, sys.nbus);  // Ybus pointer
  
  Polar** ybs = CreatePolarYbsPointer(ybs, sys.nbus);
   
  double* Ppu = createDoubleArray(Ppu, sys.nbus);   
  // PrintDoubleArray(Ppu, sys.nbus);
  
  double* Qpu = createDoubleArray(Qpu, sys.nbus);
  // PrintDoubleArray(Qpu, sys.nbus);
  
  double* Pcal = createDoubleArray(Pcal, sys.nbus);
  // PrintDoubleArray(Pcal, sys.nbus);
	
  double* Qcal = createDoubleArray(Qcal, sys.nbus);
  // PrintDoubleArray(Qcal, sys.nbus);
  
  double* DelP =  createDoubleArray(DelP, sys.nbus); 
  // PrintDoubleArray(DelP, sys.nbus);
  
  double* DelQ =  createDoubleArray(DelQ, sys.nbus);
  // PrintDoubleArray(DelQ, sys.nbus);
  
  double** J1 = createAdjMatrixN(J1, sys.nbus);
  // printAdjMatrixN(J1, sys.nbus);
  
  double** J2 = createAdjMatrixN(J2, sys.nbus);	
  // printAdjMatrixN(J2, sys.nbus);
  
  double** J3 = createAdjMatrixN(J3, sys.nbus);
  // printAdjMatrixN(J3, sys.nbus);
  
  double** J4 = createAdjMatrixN(J4, sys.nbus);	
  // printAdjMatrixN(J4, sys.nbus);
  
  Jsize = (sys.nbus - 1)*2;
    
  // Overall J Matrix
  double** J = createAdjMatrixN(J, Jsize);	// printAdjMatrixN(J, Jsize);
  
  double *x = createDoubleArray(x, Jsize);	// known vector with delp and delq combined
  // PrintDoubleArray(x, Jsize);
  
  // Sub Matrix with reduced Size to detemine the inverse of the matrix
  double** J11 = createAdjMatrixN(J11, sys.nbus-1);	// nbus -1 indicates the dimension of Jacobian 1
  // printAdjMatrixN(J11, sys.nbus-1);
  
  double** J22 = createAdjMatrixN(J22, sys.nbus-1);	// nbus -1 indicates the dimension of Jacobian 2
  // printAdjMatrixN(J22, sys.nbus-1);
  
  double** J33 = createAdjMatrixN(J33, sys.nbus-1);	// nbus -1 indicates the dimension of Jacobian 3
  // printAdjMatrixN(J33, sys.nbus-1);
  
  double** J44 = createAdjMatrixN(J44, sys.nbus-1);	// nbus -1 indicates the dimension of Jacobian 4
  // printAdjMatrixN(J44, sys.nbus-1);

  double* DelP1 =  createDoubleArray(DelP1, sys.nbus-1); 
  // PrintDoubleArray(DelP1, sys.nbus-1);
  
  double* DelQ1 =  createDoubleArray(DelQ1, sys.nbus-1); 
  // PrintDoubleArray(DelQ1, sys.nbus-1);

  double *Jele = createDoubleArray(Jele, (Jsize*Jsize));
  // PrintDoubleArray(Jele, (Jsize*Jsize));
  
  double *Sol = createDoubleArray(Sol,Jsize);
  // PrintDoubleArray(Sol, Jsize);
  
  Rect* Volt = CreateRectDoubleArray(Volt, sys.nbus);
  // printRectArray(Volt, sys.nbus);
  
  Polar* V = CreatePolarStructure(V, sys.nbus);	// Polar form of voltage
  // PrintPolarForm(V, sys.nbus);
  
  FILE *fpo1 = FilePointerCreation(fpo1, fpd.output1);
  
  FILE *fpo2 = FilePointerCreation(fpo2, fpd.output2);
  FILE *fpo3 = FilePointerCreation(fpo3, fpd.output3);
  FILE *fpo4 = FilePointerCreation(fpo4, fpd.output4);
  FILE *fpo5 = FilePointerCreation(fpo5, fpd.output5);
  FILE* fp2 = FilePointerToRead(fp2, fpd.inputfile2);
  FILE* fpo6 = FilePointerCreation(fpo6, fpd.output6);
  
  fprintf(fpo5, "*************************************************************************************\n");
  fprintf(fpo5, "Hours\t Load_Level\tPloss(kW)\tQloss(kVAr)\tBus_No\t\tVmin");
  fprintf(fpo5, "\n*************************************************************************************\n");
  
    
  Rect* Irect = CreateRectDoubleArray(Irect, sys.nbus);	// Current in Rectangular form
  // printRectArray(Irect, sys.nbus);
    
  Rect** Iflow = createRectPointer(Iflow, sys.nbus); 	// Complex Current flow pointer
  // printRectMatrix(Iflow, sys.nbus);
    
  Rect** S = createRectPointer(S, sys.nbus);  		// Complex Power flow pointer
  // printRectMatrix(S, sys.nbus);
    
  Rect* LineLoss = CreateRectDoubleArray(LineLoss, sys.nline); // Line loss pointer
  // printRectArray(LineLoss, sys.nline);	
    
  Rect* Sinj = CreateRectDoubleArray(Sinj, sys.nline);  	// Complex power injection pointer
  // printRectArray(Sinj, sys.nline);
  
  // ****** variables for Annual Load Data *****
  AnnualLoadData* ALD = CreateAnnualLoadDataPointer(ALD, hrs);
  // PrintAnnualLoadData(ALD, hrs);
  
  ALD = ReadAnnualLoadData(ALD, fp2, hrs);
  // PrintAnnualLoadData(ALD, hrs);
  
  // ************ Result Pointer Starts ************
  Result* Res = CreateResultPointer(Res, hrs);
  SystemActualLoad* ALoad = CreateSystemActLoadPointer(ALoad, one);
  
  ALoad = ComputeActualLoad(ALoad, sys, bus, one);

  // ****** Load Flow Computation starts *****
  ybus =  FormationYbus(ybus, y1, y2, line, sys); 
  // printf("\nYbus in Rectangular Form\n"); printYbusMatrix(ybus, sys.nbus); 

  ybs = ComputeMagnitudeandAngleDP(ybs, ybus, sys.nbus);    // ybus in polar form
  // printf("Ybus in polar form\n"); printYbusMatrixPolar(ybs, sys.nbus);

  for(ldl = 0; ldl < hrs; ldl++) {     // Annual load Pattern calculation
    bus = UpdatePower(bus, ALD[ldl].LoadLvl, sys);
    // printBusData(bus, sys.nbus);

    // printf("\nLoad Level no = %d\t Load Level value = %lf \n", ii+1, ALD[ii].LoadLvl);

    Ppu = CalculatepuRealPower(Ppu, bus, sys.nbus, sys.BaseMVA);   
    // printf("\nPerunit real power\n"); PrintDoubleArray(Ppu, sys.nbus);

    Qpu = CalculatepuReacPower(Qpu, bus, sys.nbus, sys.BaseMVA); 
    // printf("\nPer Unit Reactive Power\n"); PrintDoubleArray(Qpu, sys.nbus); 

    // File Printing
    fpo1 = FilePointertoWrite(fpo1, ybus, ybs, sys.nbus);

    V = LoadFlow(V, ybs, bus, Ppu, Qpu, Volt, Pcal, Qcal, J1, J2, J3, J4, DelP, DelQ, J11, J22, J33, J44, Jsize, J, x, DelP1, DelQ1, fpo3, Jele, Sol, fpo2, fpo4, sys);

    // printf("\nAfter Convergence\n"); printDoubleArray(V, sys.nbus); 

    Volt = initializeRectDoubleArray(Volt, sys.nbus);
    // printRectArray(Volt, sys.nbus);

    Volt = ConvertPoltoRect(Volt, V, sys.nbus);
    // printf("\nVoltage in Rectangular Form\n"); printRectArray(Volt, sys.nbus);    

    Irect = initializeRectDoubleArray(Irect, sys.nbus);

    Irect = ComputeBusInjection(Irect, Volt, ybus, sys.nbus);	// Bus current Injection Calculation

    // printf("\n\nCurrent in Rectangular Form\n"); printRectArray(Irect, sys.nbus);

    Iflow = initializeRectMatrix(Iflow, sys.nbus);

    Iflow = ComputeCurrentLineFlow(Iflow, line, Volt, ybus, sys.nline) ;
    // printf("\n\nCurrent Flow in Rectangular Form\n"); printRectMatrix(Iflow, sys.nbus);

    S = initializeRectMatrix(S, sys.nbus);

    S = ComputeLinePowerFlow(S, line, Iflow, Volt, ybus, sys);
    // printf("\n\nComplex Power in Rectangular Form\n"); printRectMatrix(S, sys.nbus);

    LineLoss = initializeRectDoubleArray(LineLoss, sys.nbus);

    LineLoss = CalculateLineLoss(LineLoss, S, sys, line);
    // printf("\n\nLineloss Calculated\n"); printRectArray(LineLoss, sys.nline);

    Sinj = initializeRectDoubleArray(Sinj, sys.nbus);

    Sinj = CalculateSInjection(Sinj, Volt, Irect, sys);
    // printf("\n\nComplex power injected\n");  printRectArray(Sinj, sys.nbus);

    fpo5 = FileptrLoadFlowResults(fpo5, Volt, Sinj, bus, line, S, LineLoss, V, sys, ALD, ldl, Res); 

    bus = copybusPointerData(bus, tbus, sys);   // Reset of bus pointer
    // printBusData(bus, sys.nbus);
  }
  fprintf(fpo5, "*************************************************************************************");
  
  fpo6 = SummaryFileWritten(fpo6, Res, hrs, ALoad, one, sys);	// Summary write up
  
  printf("\nTest System consider is %d bus distribution system\n", sys.nbus);
  printf("\n\nProgram Completed!!\n\nPlease check the folder to see the load flow results\n\n");
  
      
  // *** Memory free ***
  fclose(fp1);
  free(line); free(tline); free(bus); free(tbus); 
  for(i = 0; i < sys.nbus; i++) {
     free(y1[i]); free(y2[i]); } free(y1); free(y2); 
  for(i = 0; i < sys.nbus; i++) {
    free(ybus[i]); free(ybs[i]); } free(ybus); free(ybs);
  fclose(fpo1); 
  free(Ppu); free(Qpu); free(Pcal); free(Qcal); free(DelP); free(DelQ);
  for(i = 0; i < sys.nbus; i++) {
    free(J1[i]); free(J2[i]); free(J3[i]); free(J4[i]); }
  free(J1); free(J2); free(J3); free(J4);
  for(i = 0; i < Jsize; i++) {
    free(J[i]); } free(J); free(x);
  for(i = 0; i < sys.nbus-1; i++) {
    free(J11[i]); free(J22[i]); free(J33[i]); free(J44[i]); }
  free(J11); free(J22); free(J33); free(J44);
  free(DelP1); free(DelQ1); free(Jele); free(Sol); free(Volt); free(V);
  fclose(fpo2); fclose(fpo3); fclose(fpo4); fclose(fpo5);  fclose(fpo6);
  free(Irect); 
  for(i = 0; i < sys.nbus; i++) {
    free(Iflow[i]); free(S[i]); } free(Iflow); free(S);
  free(LineLoss); free(Sinj); free(ALD); fclose(fp2); free(Res);
  free(ALoad);  
  return 0;
}
