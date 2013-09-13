//////////////////////////////////////////////////////////////////////////////
//
// This program indexes diffraction spots, 
//
//
// Changes:
// V1.5: - include a fast mode: by putting large steps for position only the middle
//         of the sample is used as position. Also the plane normal is put in the 
//         middle of the 2 extremes as derived from the spotposition in relation 
//         to the ideal 2-theta ring. 
//         (in the past either 1 of the extremes of chooses for big steps) 
//       - Use box instead of eta ranges to define the spots to use (and to skip)
// v1.6: - move code from main to function
// v1.7: - split spots matrix into matrices for each ring (faster for searching)
// v1.8: - spotdata is split into omega blocks: each omega block is sorted on eta. 
//         (faster searching). Input spot data should be sorted on ringno, omegablock, eta.
//          Note: within the omegablock: spots are not sorted on omega. 
// v1.9: - clean up of code.
//       - skip planenormals if the no of matches is low (< 50% of expected)
// v2.0: - skip positions in sample if the no of mathces is low (< 50% of exp)
// v2.1: - skip only in ideal spot and pos in sample, not rotation.
// v2.3: - make (small) bins for ringnr, omega and eta. Spots are assigned to bins
//         for fast searching.
// v2.4: - removed copy of theoretical spots: matrix T. Now just add a few columns.
//       - use new input file: has build in eta and omega ranges. THis has changed comparespots as well!
// v2.5: - use a dynamic array instead of static. (also new parameters). hkls from file.
// v2.6: - make hkls file as parameter. check for eta and omega margins in comparespots (now cached)
// v2.7: - support also BCC (until now only FCC) planes
//       - unit cell structure (bcc, fcc) is now given in input file (GrainStruct) (no hkl file anymore)
//       - Change of name in input file: GrainRad -> GrainRadius    
// v2.8: - Calculate IA for the spots, as well as avg IA for the grains (output in 2 files)
// v2.9: - Include use of friedel pair: this reduces the search surface to only a few lines.
//         In case no friedel pair can be found it will use the old method: probe the whole surface.
// v2.91: - IA is added to the output files (at the last column) for indiviual spots and grains (avg)
//        - New parameter: UseFriedelPairs (1=y, 0=no)
// v3.0: - output also the best grain + spots for each starting spot (file: bestGrain_)
// v3.1: - new parameter: OutputFolder: all files are put in this folder.
// v3.2: - add output file with for EACH starting spot the score (filename: SpotScores.txt. cols: spotsid, exp, obs, fraction)
// v3.3: - change precision from float to double (IA was not calulated accurate enough, 
//         probably due to acos function). Now use type: 'RealType' everywhere
// v3.4: - solved bug: ringradii and ringnumbers in input file now consistent. Also added Boxes in params_read file.
// v3.5: - Friedelpair uses now same margins for omega and radial as margins for spots (from parameter file)
// v3.6: - OmegaRanges (and box sizes) was limited to 2, now a 72
// v3.7: - Include mixed friedelpair: reduces the amount of plane normals to search
//         if mixed friedelpair is not found: does full search on the surface.
// v3.8: - changed memory allocation slightly.
// v3.9: - changed to output only the best spot file.
// v4.0: - outputs .csv file to be used later.
//
//

// uncomment next line to get some debug info
//#define DEBUG

////////////////////////////////////////////////////////////////////////////////
// Columns of main matrices:
//
// Observed spots:
//     0. ys  [micron]     
//     1. zs  
//     2. omega    [degrees]
//     3. radius   [micron]
//     4. spot_id
//     5. ringnr
//     6. eta      [degrees]
//     7. ttheta   [degrees]
//     8. diff_idring (difference with ideal ring)  [micron]    (added in the pgm)
//
// Theoretical spots:
//     0. OrientID
//     1. spotid
//     2. indexhkl
//     3. distance  [micron]
//     4. yl        [micron]
//     5. zl        [micron]
//     6. omega     [degrees]
//     7. eta       [degrees]
//     8. theta     [degrees]
//     9. ringnr 
//     10. ys_displ [micron]
//     11. zs_displ [micron]
//     12. eta_dipl [deg]
//     13. diff_idring [micron]
//
// Output: GrainMatches
//    0-8: orientation matrix (1x9)
//    9-11: grain position xyz
//    12: No of theor spots
//    13: No of spots found
//    14: matchNo (grainNo)
//    15: avg IA
//
// Output: GrainsSpots
//    0 : spotmatch no (see GrainMatches col 14)
//    1 : dummy  
//    2 : Theor spotpos y
//    3 : Obs spotpos y
//    4 : Difference y obs-theor
//    5 : Theor spotpos z
//    6 : Obs spotpos z
//    7 : Difference z obs-theor
//    8 : Theor spotpos omega
//    9 : Obs spotpos omega
//    10 : Difference omega obs-theor
//    11 : Ref radius
//    12 : Obs radius 
//    13 : Difference radius obs-ref
//    14 : spotid
//    15 : grainmatchNo
//    16 : IA spot


#include <stdio.h>
#include <math.h> 
#include <stdlib.h> 
#include <time.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>

#define RealType double       // use 1 realtype in the whole program. 
                              // Note: single single precision turned out tp give some problems
                              // most notably with acos()     

// conversions constants
#define deg2rad 0.0174532925199433
#define rad2deg 57.2957795130823

// max array sizes    
#define MAX_N_SPOTS 2000000   // max nr of observed spots that can be stored
#define MAX_N_STEPS 1000      // Max nr of pos steps, when stepping along the diffracted ray
#define MAX_N_OR 36000        // max nr of trial orientations that can be stored (360/0.01);
#define MAX_N_MATCHES 10000   // max nr of grain matches for 1 spot
#define MAX_N_RINGS 100       // max nr of rings that can be stored (applies to the arrays ringttheta, ringhkl, etc)
#define MAX_N_HKLS 1000       // max nr of hkls that can be stored
#define MAX_N_OMEGARANGES 72  // max nr of omegaranges in input file (also max no of box sizes)

#define N_COL_THEORSPOTS 14   // number of items that is stored for each calculated spot (omega, eta, etc)
#define N_COL_OBSSPOTS 9      // number of items stored for each obs spots
#define N_COL_GRAINSPOTS 17   // nr of columns for output: y, z, omega, differences for spots of grain matches
#define N_COL_GRAINMATCHES 16 // nr of columns for output: the Matches (summary) 

// Globals
RealType ObsSpotsLab[MAX_N_SPOTS][N_COL_OBSSPOTS];              // spots, converted to lab coord
int n_spots = 0; // no of spots in global var ObsSpotsLab
 
// To store the orientation matrices
RealType OrMat[MAX_N_OR][3][3];                       

// hkls to use
int hkls[MAX_N_HKLS][4];  // columns: h k l ringno
int n_hkls = 0;

// 4d and 3d arrays for storing spots in bins. For fast lookup.
// - data[iRing][iEta][iOme] points to an array (a bin). It contains the rownumbers [0-based] of the spots in ObsSpotsLab matrix.
// - ndata holds for each bin how many spots are stored (ndata[iRing][iEta][iOme])
// - maxndata contains for each bin the capacity of the bin (maxndata[iRing][iEta][iOme]) 
// 
// ie: data[2][10][12] -> [4, 13] means spots ObsSpotsLab[4][*] and 
// ObsSpotsLab[13][*] fall into bin with: ring=3 (! not 2!), eta-segment = 10 and 
// omega-segment = 12.   
//
int ****data;
int ***ndata;
int ***maxndata;             

// the number of elements of the data arrays above
int n_ring_bins;  
int n_eta_bins;
int n_ome_bins;

// the binsizes used for the binning
RealType EtaBinSize = 0;
RealType OmeBinSize = 0;                                         
                                                       
// some macros for math calculations
#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

#define dot(v,q) \
   ((v)[0] * (q)[0] + \
    (v)[1] * (q)[1] + \
 	 (v)[2] * (q)[2])
 	  
#define CalcLength(x,y,z) sqrt((x)*(x) + (y)*(y) + (z)*(z))   	  

      
////////////////////////////////////////////////////////////////////////////////
// get the spot row numbers of a certain bin defined by ringno, eta, omega.
// eta and omega between -180 and +180 degrees.
//
int
GetBin(
  int ringno,
  RealType eta,
  RealType omega,
  int **spotRows, 
  int *nspotRows)
  
{
   // now check if this spot is in the dataset (using the bins only, a bit crude but fast way )
   int iRing, iEta, iOme, iSpot;
   iRing = ringno-1;
   iEta = floor((180+eta)/EtaBinSize);
   iOme = floor((180+omega)/OmeBinSize);
               
    int nspots = ndata[iRing][iEta][iOme];
    *spotRows = malloc(nspots*sizeof(**spotRows));
    
    if (spotRows == NULL ) {
         printf("Memory error: could not allocate memory for spotRows matrix. Memory full?\n");
         return 1;
    }
        
    // calc the diff. NOte: smallest diff in pos is choosen
    for ( iSpot = 0 ; iSpot < nspots ; iSpot++ ) {
        (*spotRows)[iSpot] = data[iRing][iEta][iOme][iSpot];
    }
            
    *nspotRows = nspots;
    return 0;    
}


  
////////////////////////////////////////////////////////////////////////////////  
// Finds a value in a column of a matrix 
// returns the rowno. -1 if the value is not found.
// 
// example: FindInMatrix(&Mat[0][0], nrows, ncols, 2, 10.1, &SpotRowNo);
//    searches for 10.1 in column-index 2 (0 based, so actually column 3!)
//
void  
FindInMatrix(
  RealType *aMatrixp,
  int nrows,
  int ncols,
  int SearchColumn,  
  RealType aVal,
  int *idx)
   
{
  int r, LIndex;
  *idx = -1;
  
  for (r=0 ; r< nrows ; r++) {
    LIndex = (r*ncols) + SearchColumn;
    if (aMatrixp[LIndex] == aVal ) {
       *idx = r;
       break;
    }
  }
}



void PrintSpotID(RealType SpotID)
{
    int SpotRowNo;
    
      // find row for given spotID 
      FindInMatrix(&ObsSpotsLab[0][0], n_spots, N_COL_OBSSPOTS, 4, SpotID, &SpotRowNo);
    
      if (SpotRowNo == -1) {
         return;
      }
      
      RealType ys     = ObsSpotsLab[SpotRowNo][0];   
      RealType zs     = ObsSpotsLab[SpotRowNo][1];
      RealType omega  = ObsSpotsLab[SpotRowNo][2];  
      //RealType Rad = ObsSpotsLab[SpotRowNo][3];
      RealType eta    = ObsSpotsLab[SpotRowNo][6];
      RealType ttheta = ObsSpotsLab[SpotRowNo][7];
      int   ringnr = (int) ObsSpotsLab[SpotRowNo][5];  
      
      printf("%8s %10s %9s %9s %9s %9s %9s %7s\n", "SpotID", "SpotRowNo", "ys", "zs", "omega", "eta", "ttheta", "ringno");
      printf("%8.0f %10d %9.2f %9.2f %9.3f %9.3f %9.3f %7d\n\n", SpotID, SpotRowNo, ys, zs, omega, eta, ttheta, ringnr);
}


////////////////////////////////////////////////////////////////////////////////
// allocates 2d array
// returns NULL if failed.
RealType** 
allocMatrix(int nrows, int ncols)
{
   RealType** arr;
   int i;

   arr = malloc(nrows * sizeof(*arr));
   if (arr == NULL ) {
      return NULL;
   }      
   for ( i = 0 ; i < nrows ; i++) {
      arr[i] = malloc(ncols * sizeof(*arr[i]));
      if (arr[i] == NULL ) {
         return NULL;
     }     
   }
   
   return arr;
}


////////////////////////////////////////////////////////////////////////////////
void
FreeMemMatrix(
   RealType **mat,
   int nrows)
   
{
   int r;
   
   for ( r = 0 ; r < nrows ; r++) {
      free(mat[r]);
   }
   
   free(mat);
}


////////////////////////////////////////////////////////////////////////////////
void
InitArrayI(
  int anArray[],
  int nel)
  
{
   memset(anArray, 0, sizeof(int)*nel);   

//   int i;  
//    for (i=0 ; i < nel ; ++i)
//    {
//      anArray[i] = 0;
//    } 
}


////////////////////////////////////////////////////////////////////////////////
// NB: only works for matrices that are contiguous!
//
void
InitMatrixCI(
  int * aMatrix,
  int nrows,
  int ncols)
  
{
   memset(aMatrix, 0, sizeof(int)*nrows*ncols);

//    int r,c;
//    for (r=0 ; r < nrows ; ++r) {
//       for (c=0 ; c < ncols ; ++c) {
//          aMatrix[r * ncols + c] = 0;
//       }
//    } 
}
    
    
////////////////////////////////////////////////////////////////////////////////
// NB: only works for matrices that are contiguous!
void
InitMatrixCF(
   RealType *aMatrix,
   int nrows,
   int ncols)
  
{
   memset(aMatrix, 0, sizeof(RealType)*nrows*ncols);
   
//    int r,c;
//   
//    for (r=0 ; r < nrows ; ++r) {
//       for (c=0 ; c < ncols ; ++c) {
//          aMatrix[r * ncols + c] = 0;
//       }
//    } 
}
    
    
////////////////////////////////////////////////////////////////////////////////
void
PrintMatrix3(
   int aMatrix[][3],
   int nrows) {
    
   int r, c;
   for ( r = 0; r<nrows ; r++) {
      for ( c = 0; c<3 ; c++) {
         printf("%d ", aMatrix[r][c]);
      }
      printf("\n");
   }
}


////////////////////////////////////////////////////////////////////////////////
// writes a matrix with ints to file.
//
int
WriteMatrixIp(
   char FileName[],
   int **aMatrixp,
   int nrows,
   int ncols) 
  
{
   int r, c;
   FILE *fp;
  
   fp = fopen(FileName, "w");
    
   if (fp==NULL) {
      printf("Cannot open file, %s\n", FileName);
      return (1);
   }
    
   for(r=0; r<nrows; r++) { 
      for(c=0;c<ncols; c++) {
         fprintf(fp, "%d ", aMatrixp[r][c]);
      }
      fprintf(fp, "\n");
   }
  
   fclose(fp);
   return(0);
}



////////////////////////////////////////////////////////////////////////////////
// writes matrix with RealType values, with optional header
//
void
WriteMatrixWithHeaderFp(
   char FileName[],
   RealType **aMatrixp,
   int nrows,
   int ncols,
   char header[])  // give empty string ("\0") to skip header) 

{
   int r, c;
   FILE *fp;
  
   fp = fopen(FileName, "w");
    
   if (fp==NULL) {
      printf("Cannot open file, %s\n", FileName);
      return;
   }
  
   // write header
   if (header[0] != '\0')   {
      fprintf(fp, "%s\n", header);
   }
  
   // write data
   for(r=0; r<nrows; r++) { 
      for(c=0;c<ncols; c++) {
         fprintf(fp, "%14f ", aMatrixp[r][c]);         
      }
      fprintf(fp, "\n");
   }
  
   fclose(fp);
}

////////////////////////////////////////////////////////////////////////////////
// writes a matrix with RealTypes to file.
//
void
WriteMatrixFp(
   char FileName[],
   RealType **aMatrixp,
   int nrows,
   int ncols) 
  
{
   WriteMatrixWithHeaderFp(FileName, aMatrixp, nrows,  ncols, "");
}



////////////////////////////////////////////////////////////////////////////////
//  write output file with the spots found for a grain.
//
int
WriteMatrixGrainsSpots(
   char FileName[],
   RealType **aMatrix,
   int nrows,
   int ncols)

{
   int r, c;
   FILE *fp;
  
   fp = fopen(FileName, "w");
    
   if (fp == NULL) {
      printf("Cannot open file, %s\n", FileName);
      return 1;
   }
  
   // write header
   char header[] = "No Dummy YCalc YObs YDiff ZCalc ZObs ZDiff WCalc WObs WDiff RadRef RadObs RadDiff SpotId MatchNo IA";
   fprintf(fp, "%s\n", header);
  
   // write data
   for(r=0; r<nrows; r++) { 
      for(c=0; c<ncols; c++) {
         fprintf(fp, "%lf ", aMatrix[r][c]);
      }
      fprintf(fp, "\n");
   }
  
   fclose(fp);
   return 0;   
}


////////////////////////////////////////////////////////////////////////////////
// NB only works for 3d matrices that are contiguous in memory!
void
WriteMatrix3DFp(
  char FileName[],
  RealType *aMatrixp,
  int nlayers, 
  int nrows,
  int ncols) {

  int l, r, c;
  FILE *fp;
  
  fp = fopen(FileName, "w");
    
  if (fp==NULL) {
    printf("Cannot open file, %s\n", FileName);
    return;
  }
    
  for (l=0; l<nlayers; l++) {
    for(r=0; r<nrows; r++) { 
      for(c=0; c<ncols; c++) {
        fprintf(fp, "%lf ",aMatrixp[l*ncols*nrows + r * ncols + c]);
      }
      fprintf(fp, "\n");
    }
  }
  
  fclose(fp);
}



////////////////////////////////////////////////////////////////////////////////
void
WriteArrayF(
   char FileName[],
   RealType anArray[],
   int nel)  {
  
   int i;
   FILE *fp;
  
   fp = fopen(FileName, "w");
    
   if (fp==NULL) {
      printf("Cannot open file: %s\n", FileName);
      return;
   }   
  
   for (i = 0; i<nel ; i++) {
      fprintf(fp, "%lf\n", anArray[i]);
   }
  
   fclose(fp);  
}


////////////////////////////////////////////////////////////////////////////////
void
WriteArrayI(
   char FileName[],
   int anArray[],
   int nel)  {
  
   int i;
   FILE *fp;
  
   fp = fopen(FileName, "w");
    
   if (fp==NULL) {
      printf("Cannot open file: %s\n", FileName);
      return;
   }  
  
   for (i = 0; i<nel ; i++) {
      fprintf(fp, "%d\n", anArray[i]);
   }
  
   fclose(fp);  
}


////////////////////////////////////////////////////////////////////////////////
void
PrintMatrixFp(
   RealType **aMatrixp,
   int nrows,
   int ncols,
   int linenr) {
    
   int r, c;
   for ( r = 0; r<nrows ; r++) {
      if (linenr == 1) printf("%d ", r+1);
      for ( c = 0; c<ncols ; c++) {
         printf("%lf ", aMatrixp[r][c]);
      }
      printf("\n");
   }
}


////////////////////////////////////////////////////////////////////////////////
// for printing a 2d matrix, treating it like a long 1d matrix. 
// so the first argument is a pointer
// call with PrintMatrix1Df(&matrix[0][0]...)
// NB: only works for matrices that are contiguous in memory.
void
PrintMatrixCF(
   RealType *aMatrixp,
   int nrows,
   int ncols,
   int linenr) {
    
   int r, c;
   for ( r = 0; r<nrows ; r++) {
      if (linenr == 1) printf("%d ", r+1);
      for ( c = 0; c<ncols ; c++) {
         printf("%lf ", aMatrixp[r * ncols + c]); 
      }
      printf("\n");
   }
}



////////////////////////////////////////////////////////////////////////////////
void
PrintMatrixF(
   RealType aMatrix[][3],
   int nrows) {
    
   int r, c;
   for ( r = 0; r<nrows ; r++) {
      printf("%d ", r);
      for ( c = 0; c<3 ; c++) {
         printf("%lf ", aMatrix[r][c]);
      }
      printf("\n");
   }
}


////////////////////////////////////////////////////////////////////////////////
void
PrintArray(
   int anArray[],
   int nel)  
{
   int i;
  
   for (i = 0; i<nel ; i++) {
      printf("%d ", anArray[i]);
   }
   printf("\n");
}


////////////////////////////////////////////////////////////////////////////////
void
PrintArrayF(
   RealType anArray[],
   int nel)  
{
   int i;
  
   for (i = 0; i<nel ; i++) {
      printf("%lf ", anArray[i]);
   } 
   printf("\n");
}


////////////////////////////////////////////////////////////////////////////////
void
PrintArrayVertF(
   RealType anArray[],
   int nel)  {
  
   int i;
  
   for (i = 0; i<nel ; i++) {
      printf("%lf\n", anArray[i]);
   }
}


////////////////////////////////////////////////////////////////////////////////
void
MatrixMultF33(
   RealType m[3][3],
   RealType n[3][3],
   RealType res[3][3]) 

{
   int r;
   
   for (r=0; r<3; r++) {
      res[r][0] = m[r][0]*n[0][0] +
                  m[r][1]*n[1][0] +
                  m[r][2]*n[2][0];

      res[r][1] = m[r][0]*n[0][1] +
                  m[r][1]*n[1][1] +
                  m[r][2]*n[2][1];
                   
      res[r][2] = m[r][0]*n[0][2] +
                  m[r][1]*n[1][2] +
                  m[r][2]*n[2][2];
                    
   }    
}


////////////////////////////////////////////////////////////////////////////////
void
MatrixMultF(
   RealType m[3][3],
   RealType v[3],
   RealType r[3])  

{
   int i;
   
   for (i=0; i<3; i++) {
      r[i] = m[i][0]*v[0] +
             m[i][1]*v[1] +
             m[i][2]*v[2]; 
                    
   }    
} 


////////////////////////////////////////////////////////////////////////////////
void
MatrixMult(
   RealType m[3][3],
   int  v[3],
   RealType r[3]) 

{

   int i;
   
   for (i=0; i<3; i++) {
      r[i] = m[i][0]*v[0] +
             m[i][1]*v[1] +
             m[i][2]*v[2]; 
   }    
} 


////////////////////////////////////////////////////////////////////////////////
RealType
min(RealType a, RealType b)
{
    return (a < b ? a : b);
}


////////////////////////////////////////////////////////////////////////////////
RealType
max(RealType a, RealType b)
{
    return (a > b ? a : b);
}

      
      	

////////////////////////////////////////////////////////////////////////////////
int
ReadHKLs( char FileName[])

{
   FILE *fp;
   n_hkls = 0;   
   int h, k, l, ringno;
   
   fp = fopen(FileName, "r"); // 
   if ( fp == NULL ) {
      printf("Error: cannot open file: %s.\n", FileName);
      return (1);
   }
   
   while ( fscanf(fp, "%d %d %d %d", &h, &k, &l, &ringno) != EOF ) {
   //while ( fscanf(fp, "%d %d %d", &h, &k, &l) != EOF ) {   
      if (n_hkls >= MAX_N_HKLS)  {
         printf("Error: number of hkl's in file %s exceeds maximum (%d).\n", FileName, MAX_N_HKLS);
         return 2;
         break;
      }   
      hkls[n_hkls][0] = h;
      hkls[n_hkls][1] = k;
      hkls[n_hkls][2] = l;      
      hkls[n_hkls][3] = ringno;      
      
      n_hkls++;
   }
   fclose(fp);
   
   return (0);
}


////////////////////////////////////////////////////////////////////////////////
// returns the all hkls for given ringnumbers in FCC.
void
GenerateHKLsFCC(
  int ringnos[],
  int nrings,  
  int hkls[][4],
  int *nhkls)
  
{

#define N_HKLSFCC 58  // should be smaller then MAX_N_HKLS
int hklsFCC[N_HKLSFCC][4] ={
    {-1,-1,-1,1},
    { 1, 1, 1,1},
    {-1,-1, 1,1},
    {-1, 1,-1,1},
    {-1, 1, 1,1},
    { 1,-1,-1,1},
    { 1,-1, 1,1},
    { 1, 1,-1,1},
    
    {-2, 0, 0,2},  
    { 0,-2, 0,2},
    { 0, 0,-2,2},  
    { 0, 0, 2,2},  
    { 0, 2, 0,2},  
    { 2, 0, 0,2},
      
    {-2,-2, 0,3},
    {-2, 0,-2,3},
    {-2, 0, 2,3},
    {-2, 2, 0,3},
    { 0,-2,-2,3},
    { 0,-2, 2,3},
    { 0, 2,-2,3},
    { 0, 2, 2,3},
    { 2,-2, 0,3},
    { 2, 0,-2,3},
    { 2, 0, 2,3},
    { 2, 2, 0,3},
    
    {-3,-1,-1,4},
    {-1,-3,-1,4},
    { 1, 3, 1,4},
    { 3, 1, 1,4},
    {-3,-1, 1,4},
    {-3, 1,-1,4},
    {-3, 1, 1,4},
    {-1,-3, 1,4},
    {-1,-1,-3,4},
    {-1,-1, 3,4},
    {-1, 1,-3,4},
    {-1, 1, 3,4},
    {-1, 3,-1,4},
    {-1, 3, 1,4},
    { 1,-3,-1,4},
    { 1,-3, 1,4},
    { 1,-1,-3,4},
    { 1,-1, 3,4},
    { 1, 1,-3,4},
    { 1, 1, 3,4},
    { 1, 3,-1,4},
    { 3,-1,-1,4},
    { 3,-1, 1,4},
    { 3, 1,-1,4},
    
    {-2,-2,-2,5},
    { 2, 2, 2,5},
    {-2,-2, 2,5},
    {-2, 2,-2,5},
    {-2, 2, 2,5},
    { 2,-2,-2,5},
    { 2,-2, 2,5},
    { 2, 2,-2,5} };
 
   int i, ringno, rowno;
   
   *nhkls = 0;
   // copy only the rings needed:  check all the rings
   for ( i = 0 ; i < nrings ; i++ ) {
      ringno = ringnos[i];
      
      // check all hkls in the database
      for ( rowno = 0 ; rowno < N_HKLSFCC ; rowno++ ) {
          if (hklsFCC[rowno][3] == ringno) { 
             hkls[*nhkls][0] = hklsFCC[rowno][0]; 
             hkls[*nhkls][1] = hklsFCC[rowno][1]; 
             hkls[*nhkls][2] = hklsFCC[rowno][2];
             hkls[*nhkls][3] = hklsFCC[rowno][3];
             (*nhkls)++;             
          }
      }
   }
}

    
////////////////////////////////////////////////////////////////////////////////
// generate for given ringno's all the hkls for bcc
//
void
GenerateHKLsBCC(
  int ringnos[],
  int nrings,  
  int hkls[][4],
  int *nhkls)
  
{

#define N_HKLSBCC 78    // should be smaller then MAX_N_HKLS
int hklsBCC[N_HKLSBCC][4] ={
    {-1, -1,  0, 1}, 
    {-1,  0, -1, 1},
    {-1,  0,  1, 1},
    {-1,  1,  0, 1},
    { 0, -1, -1, 1},
    { 0, -1,  1, 1},
    { 0,  1, -1, 1},
    { 0,  1,  1, 1},
    { 1, -1,  0, 1},
    { 1,  0, -1, 1},
    { 1,  0,  1, 1},
    { 1,  1,  0, 1},
    {-2,  0,  0, 2}, 
    { 0, -2,  0, 2},
    { 0,  0, -2, 2}, 
    { 0,  0,  2, 2},  
    { 0,  2,  0, 2},  
    { 2,  0,  0, 2},  
    {-2, -1, -1, 3},
    {-1, -2, -1, 3},
    { 1,  2,  1, 3},
    { 2,  1,  1, 3},
    {-2, -1,  1, 3},
    {-2,  1, -1, 3},
    {-2,  1,  1, 3},
    {-1, -2,  1, 3},
    {-1, -1, -2, 3},
    {-1, -1,  2, 3},
    {-1,  1, -2, 3},
    {-1,  1,  2, 3},
    {-1,  2, -1, 3},
    {-1,  2,  1, 3},
    { 1, -2, -1, 3},
    { 1, -2,  1, 3},
    { 1, -1, -2, 3},
    { 1, -1,  2, 3},
    { 1,  1, -2, 3},
    { 1,  1,  2, 3},
    { 1,  2, -1, 3},
    { 2, -1, -1, 3},
    { 2, -1,  1, 3},
    { 2,  1, -1, 3}, 
    {-2, -2,  0, 4}, 
    {-2,  0, -2, 4},
    {-2,  0,  2, 4},
    {-2,  2,  0, 4},
    { 0, -2, -2, 4},
    { 0, -2,  2, 4},
    { 0,  2, -2, 4},
    { 0,  2,  2, 4},
    { 2, -2,  0, 4},
    { 2,  0, -2, 4},
    { 2,  0,  2, 4},
    { 2,  2,  0, 4},
    { 3,  1,  0, 5},
    { 3,  0,  1, 5},
    { 1,  3,  0, 5},    
    { 1,  0,  3, 5},
    { 0,  3,  1, 5},        
    { 0,  1,  3, 5},    
    {-3,  1,  0, 5},
    {-3,  0,  1, 5},
    { 1, -3,  0, 5},    
    { 1,  0, -3, 5},
    { 0, -3,  1, 5},        
    { 0,  1, -3, 5},    
    { 3, -1,  0, 5},
    { 3,  0, -1, 5},
    {-1,  3,  0, 5},    
    {-1,  0,  3, 5},
    { 0,  3, -1, 5},        
    { 0, -1,  3, 5},    
    {-3, -1,  0, 5},
    {-3,  0, -1, 5},
    {-1, -3,  0, 5},    
    {-1,  0, -3, 5},
    { 0, -3, -1, 5},        
    { 0, -1, -3, 5} };          
 
   int i, ringno, rowno;
   
   *nhkls = 0;
   // copy only the rings needed: check all the rings
   for ( i = 0 ; i < nrings ; i++ ) {
      ringno = ringnos[i];
      
      // check all hkls in the database
      for ( rowno = 0 ; rowno < N_HKLSBCC ; rowno++ ) {
          if (hklsBCC[rowno][3] == ringno) { 
             hkls[*nhkls][0] = hklsBCC[rowno][0]; 
             hkls[*nhkls][1] = hklsBCC[rowno][1]; 
             hkls[*nhkls][2] = hklsBCC[rowno][2];
             hkls[*nhkls][3] = hklsBCC[rowno][3];
             (*nhkls)++;             
          }
      }
   }
}



////////////////////////////////////////////////////////////////////////////////
void
CalcInternalAngle(
   RealType x1, 
   RealType y1, 
   RealType z1, 
   RealType x2, 
   RealType y2, 
   RealType z2, 
   RealType *ia)
   
{
   RealType v1[3];
   RealType v2[3];
   
   v1[0] = x1;
   v1[1] = y1;   
   v1[2] = z1;
   
   v2[0] = x2;
   v2[1] = y2;   
   v2[2] = z2;
    
   RealType l1 = CalcLength(x1, y1 ,z1);
   RealType l2 = CalcLength(x2, y2, z2);
   RealType tmp = dot(v1, v2)/(l1*l2);
   
   if (tmp > 1 ) { tmp = 1;  }  
   if (tmp < -1 ) {tmp = -1; }   
   
   *ia = rad2deg * acos(tmp);
}


////////////////////////////////////////////////////////////////////////////////
int
ReadMatrix(
   char  FileName[],
   RealType aMatrix[][30],
   int * nrows)
{
   int i, j;
   FILE *fp;
   char buffer[1000];
   
   *nrows = 0;   
   
   fp = fopen(FileName, "r"); // 
   if (fp==NULL) {
     printf("Cannot open file: %s.\n", FileName);
     return 1;
   }
  
   // skip header
   fgets(buffer, 1001, fp);   
   int EofReached = 0;
   for(i=0; i<MAX_N_SPOTS; i++) {
      for(j=0; j<30; j++) {
         if (fscanf(fp, "%lf", &aMatrix[i][j]) == EOF)     EofReached = 1; 
        
       }
       if (EofReached == 1) break;
       *nrows = *nrows + 1;            
   }
    
   fclose(fp);
   
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
//
// Read corrected spot file.
//
// See for columns of cor file top of the source file.   
//
int
ReadMatrixCor(
   char  FileName[],
   RealType aMatrix[][N_COL_OBSSPOTS],     
   int * nrows)

{
   FILE *fp;
   int maxLineLength = 10000;   
   char aline[maxLineLength];
   RealType  ys, zs, omega, radius, spotid, ringnr, eta, ttheta;
   *nrows = 0;   
   
   fp = fopen(FileName, "rt"); // 
   if (fp==NULL) {
      printf("Cannot open file: %s.\n", FileName);
      return 1;
   }
   // skip header
   fgets(aline, maxLineLength, fp);   

   // get line by line
   while(fgets(aline, maxLineLength, fp) != NULL)
   {
      // check if not to much data
      if ( *nrows >= MAX_N_SPOTS )
      {
         printf("Error: Too many spots in datafile %s (Max = %d)\n", FileName, MAX_N_SPOTS );
         fclose(fp);
         return(2);
      }
   
	   // convert the string to NoCol RealTypes
      sscanf(aline, "%lf %lf %lf %lf %lf %lf %lf %lf\
                     %*f %*f %*f %*f %*f %*f %*f",
                     &ys, &zs, &omega, &radius, &spotid, &ringnr, &eta, &ttheta);
                     
      aMatrix[*nrows][0] = ys;
      aMatrix[*nrows][1] = zs;
      aMatrix[*nrows][2] = omega;
      aMatrix[*nrows][3] = radius;
      aMatrix[*nrows][4] = spotid;
      aMatrix[*nrows][5] = ringnr;
      aMatrix[*nrows][6] = eta; 
      aMatrix[*nrows][7] = ttheta;
      (*nrows)++;                     
   }
    
   fclose(fp);
     
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
int 
ReadArrayI(
   char  FileName[],
   int aArray[],
   int * nel)
   
{
   int i;
   FILE *fp;
   
   *nel = 0;   
   
   fp = fopen(FileName, "r"); // 
   if (fp==NULL) {
      printf("Cannot open file: %s.\n", FileName);
      return (1);
   }
  
   // skip header
//   fgets(buffer, 1001, fp);   
   for(i=0; i<MAX_N_SPOTS; i++) {
      if (fscanf(fp, "%d", &aArray[i]) == EOF)  break;
      *nel = *nel + 1;       
   }
   
   fclose(fp);
   
   return (0);
}


////////////////////////////////////////////////////////////////////////////////
// 
// % Rotates a 3d vector around z over an angle alpha.
// % For right hand system (as in fable) this means +alpha = ccw.
// %
// % Input:
// %     v1:  xyz vector 
// %     alpha: angle [degrees]
// %
// % Output:
// %     rotated vector
// %
// % Date: 18-08-2010
// %
void
RotateAroundZ(
   RealType v1[3],        
   RealType alpha, 
   RealType v2[3]) 
{
   RealType cosa = cos(alpha*deg2rad);
   RealType sina = sin(alpha*deg2rad);
   
   RealType mat[3][3] = {{ cosa, -sina, 0 },
          { sina,  cosa, 0 },
          { 0, 0, 1}};
      
   MatrixMultF(mat, v1, v2);   
}




////////////////////////////////////////////////////////////////////////////////
// % Calculates eta angle of a coordinate y,z, with 0 degrees = the z axis pointing up.
// % 
// % input:
// %   y, z: coordinates with y to the right, and z up (right hand system)
// %
// % output:
// %   alpha: angle between 0 (up) and -180 ccw or 180 cw [degrees]
// %
// % Date: 18-08-2010
// % 
void
CalcEtaAngle(
   RealType y,
   RealType z, 
   RealType *alpha) {
  
   // calc angle using dot product.
   *alpha = rad2deg * acos(z/sqrt(y*y+z*z));  // rad -> deg
    
   // alpha is now [0, 180], for positive y it should be between [0, -180]
   if (y > 0)    *alpha = -*alpha;
}



////////////////////////////////////////////////////////////////////////////////
//       
// % Calculates the position of a spot on a detector in lab coordinates. 
// % The origin is the center of the sample rotation.
// %
// % Input:
// %    RingRadius : the radius of the ring of the spot [any unit, output is in same units]
// %    eta : angle along the ring [degrees] 0 degrees = up, clockwise is
// %          positive (seen from the beam).
// %
// % Output:
// %    yl,zl: spot position in lab coordinates. [unit of distance input]
// %        (x= direction of beam, y=to the left, z=up)
// %
// % Date: 23-08-2010
// %
void
CalcSpotPosition(
   RealType RingRadius, 
   RealType eta, 
   RealType *yl,
   RealType *zl)
  
{
   RealType etaRad = deg2rad * eta;
    
   *yl = -(sin(etaRad)*RingRadius); // y inversed: to the 'door' (=left) is +!
   *zl =   cos(etaRad)*RingRadius;
}       
       

////////////////////////////////////////////////////////////////////////////////
// % This function finds the rotation angles (omegas) for which diffration occurs
// %
// % Input:
// %   x y z: the coordinates of the gvector (for omega = 0)
// %          +x forward
// %          +y to the left (!)
// %          +z up
// %   theta: diffraction angle [degrees]
// % 
// % Output: 
// %   The rotation angles w for which diffraction occurs [degrees: 0-360].
// %   Positive omega means: going to the quadrant where x and y are
// %   positive. 
// %   In fable this means ccw!
// %
// % NB: for very small values of y (small but > 1e6) the routine might give
// % 4 answers (instead of 2) which lie very close together.
// %
// % Date: 18-08-2010
// %
// %

void
CalcOmega(
   RealType x,
   RealType y,
   RealType z,
   RealType theta, 
   RealType omegas[4],
   RealType etas[4],
   int * nsol) {

// % solve the equation -x cos w + y sin w = sin(theta) * len
// % it simply comes from calculating the angle between incoming beam
// % (represented by -1 0 0, note the minus sign!) and the gvec.
// %

   *nsol = 0;
   RealType ome;
   RealType len= sqrt(x*x + y*y + z*z);
   RealType v=sin(theta*deg2rad)*len;
  
   // % in case y is 0 use short version: -x cos w = v
   // % NB use radians in this part of the code: its faster
   RealType almostzero = 1e-4;
   if ( fabs(y) < almostzero ) {
      if (x != 0) {
         RealType cosome1 = -v/x;
         if (fabs(cosome1 <= 1)) {
            ome = acos(cosome1)*rad2deg;
            omegas[*nsol] = ome;
            *nsol = *nsol + 1;
            omegas[*nsol] = -ome;  // todo: not in range[0 360]
            *nsol = *nsol + 1;
         }
      }
   }    
   else { //% y != 0
      RealType y2 = y*y;
      RealType a = 1 + ((x*x) / y2);
      RealType b = (2*v*x) / y2;
      RealType c = ((v*v) / y2) - 1;
      RealType discr = b*b - 4*a*c;    
      
      RealType ome1a;
      RealType ome1b;
      RealType ome2a;
      RealType ome2b;
      RealType cosome1;
      RealType cosome2;      
      
      RealType eqa, eqb, diffa, diffb;
  
      if (discr >= 0) {
         cosome1 = (-b + sqrt(discr))/(2*a);
         if (fabs(cosome1) <= 1) {       
            ome1a = acos(cosome1);
            ome1b = -ome1a;
       
            // not all omegas found are valid! Accept answer only if it equals v!
            eqa = -x*cos(ome1a) + y*sin(ome1a);
            diffa = fabs(eqa - v);
            eqb = -x*cos(ome1b) + y*sin(ome1b);
            diffb = fabs(eqb - v);
            
            // take the closest answer
            if (diffa < diffb ) {
               omegas[*nsol] = ome1a*rad2deg;
               *nsol = *nsol + 1;
            }               
            else {
               omegas[*nsol] = ome1b*rad2deg;
               *nsol = *nsol + 1;             
            }
         }
  
         cosome2 = (-b - sqrt(discr))/(2*a);
         if (fabs(cosome2) <= 1) {    
            ome2a = acos(cosome2);
            ome2b = -ome2a;
          
            eqa = -x*cos(ome2a) + y*sin(ome2a);
            diffa = fabs(eqa - v);
            eqb = -x*cos(ome2b) + y*sin(ome2b);
            diffb = fabs(eqb - v); 
            
            if (diffa < diffb) {
               omegas[*nsol] = ome2a*rad2deg;
               *nsol = *nsol + 1;
            }
            else {
               omegas[*nsol] = ome2b*rad2deg;
               *nsol = *nsol + 1;                           
            }
         }   
      }
   }      

   // find corresponding eta's
   RealType gw[3];
   RealType gv[3]={x,y,z};
   RealType eta;
   int indexOme;
   for (indexOme = 0; indexOme < *nsol; indexOme++) {
      // first get the rotated vector, then the y and z coordinates gives
      // directly the eta angle [-180, 180]. 
      // y and z are in the detector plane
      RotateAroundZ(gv, omegas[indexOme], gw);
      CalcEtaAngle(gw[1],gw[2], &eta);
      etas[indexOme] = eta; 
   }
}


////////////////////////////////////////////////////////////////////////////////
// % This function calculates the diffraction spots for transmission
// % diffraction, while the sample rotates 360 degrees around the z-axis
// % (vertical axis).
// %
// % Input: 
// %   orientation (rotation matrix 3x3)
// %   hkl planes
// %   wavelength of x-ray beam [Angstrom]
// %   distance sample-detector [micron]
// %
// % Output:
// %   spots: For each orientation the diffraction spots (see below)
// %   spotsnr: The number of spots
// %
// %   For each spot (in this order, columnwise):
// %     OrientID   : a number from 1..number of orientations. (obsolete) 
// %     spotid     : a number from 1..number of spots (total)
// %     indexhkl   : a number that indicates the hkl plane. To get the hkl plane
// %                  use hkls(indexhkl)
// %     xl, yl, zl : the lab coordinates of the spot on the detector [micron] (x=direction of the
// %                  beam,  y=to the left, z=up, (0,0,0) is the center of rotation)
// %     omega      : rotation angle, from top seen: ccw is positive [degrees]
// %     eta        : position of spot on the ring, seen from the beam: cw is postive [deg]
// %     theta      : diffraction angle [degrees].
//       ringnr     : the ringnr of the spot
// %  
// % Date: 03-09-2010
// %
// % version _Furnace:
// %   - added 3rd and 4th ring, but only certain etas for certain omegas (special case
// %       when the furnace is used: the top and bottom are missing on 1 side)

////////////////////////////////////////////////////////////////////////////////
void   
CalcDiffrSpots_Furnace(
   RealType OrientMatrix[3][3], 
   RealType LatticeConstant, 
   RealType Wavelength , 
   RealType distance,
   RealType RingRadii[],  // the observed ring radii
   RealType OmegaRange[][2],  // omegaranges min - max [-180 .. 180]]  
   RealType BoxSizes[][4],     // For each omegerange: size of the box on the detector: all spots outside are removed.
   int NOmegaRanges,  // no of OmegaRanges and BoxSizes    
   RealType ExcludePoleAngle,
   RealType **spots,  // TODO check spots returned correctly?
   int   *nspots) 
  
{
  // get number of planes
//  printf("number of planes: %d\n", n_hkls);
  
   // calculate dspacing and theta for the planes in hkls (store in the table)
   int i, OmegaRangeNo;
   RealType DSpacings[MAX_N_HKLS];
   RealType thetas[MAX_N_HKLS];
   RealType ds;
   RealType theta; 
   int KeepSpot;
   
   // calculate theta and dspacing
   for (i = 0 ;i < n_hkls; i++) { 
      int sumhkl = hkls[i][0]*hkls[i][0] + hkls[i][1]*hkls[i][1] + hkls[i][2]*hkls[i][2];
      ds = sqrt(LatticeConstant*LatticeConstant / sumhkl);
      theta = rad2deg * asin(Wavelength/(2*ds));    
      DSpacings[i] = ds;
      thetas[i] = theta;
   }
   
   // Each spots has 1 row with information about eta, omega, etc..
   // each hkl plane gives 2 spots for 360 degrees.
   
   // use rotation matrix for given orientation
   // generate for this orientation for each hkl the spots
   int Ghkl[3]; // vector
   int indexhkl;
   RealType Gc[3]; 
   RealType omegas[4];
   RealType etas[4];
   RealType yl;
   RealType zl;
   int nspotsPlane;
   int spotnr = 0;  // spot id's for this orientation
   int spotid = 0;
   int OrientID = 0;
   int ringnr = 0;
   
       
   for (indexhkl=0; indexhkl < n_hkls ; indexhkl++)  {
      // g vector for non rotated crystal (in ref system)
      Ghkl[0] = hkls[indexhkl][0];
      Ghkl[1] = hkls[indexhkl][1];    
      Ghkl[2] = hkls[indexhkl][2];
      
      // get ringnr   
      ringnr = hkls[indexhkl][3];
      RealType RingRadius = RingRadii[ringnr];
      
      // calculate gvector for given orientation: Gc (in lab coordinates!)
      MatrixMult(OrientMatrix,Ghkl, Gc);
      
      // calculate omega angles for which diffraction occurs. Eta is
      // also calculated.
      ds    = DSpacings[indexhkl];
      theta = thetas[indexhkl];
      CalcOmega(Gc[0], Gc[1], Gc[2], theta, omegas, etas, &nspotsPlane);
      
      // calculate spot on detector in lab coordinates (xl,yl,zl)
      for (i=0 ; i<nspotsPlane ; i++) {
         RealType Omega = omegas[i];
         RealType Eta = etas[i]; 
         RealType EtaAbs =  fabs(Eta);
      
         // remove spots on the poles
         if ((EtaAbs < ExcludePoleAngle ) || ((180-EtaAbs) < ExcludePoleAngle)) continue; 
         
         CalcSpotPosition(RingRadius, etas[i], &(yl), &(zl));
         
         for (OmegaRangeNo = 0 ; OmegaRangeNo < NOmegaRanges ; OmegaRangeNo++ ) {
            KeepSpot = 0;
            // if spot inside Omegarange and inside box, keep it
            if ( (Omega > OmegaRange[OmegaRangeNo][0]) && 
                 (Omega < OmegaRange[OmegaRangeNo][1]) &&
                 (yl > BoxSizes[OmegaRangeNo][0]) && 
                 (yl < BoxSizes[OmegaRangeNo][1]) && 
                 (zl > BoxSizes[OmegaRangeNo][2]) && 
                 (zl < BoxSizes[OmegaRangeNo][3]) ) {
               KeepSpot = 1;
               break;
            }
         }
      
         if (KeepSpot) {
            spots[spotnr][0] = OrientID;
            spots[spotnr][1] = spotid;
            spots[spotnr][2] = indexhkl;
            spots[spotnr][3] = distance;  // xl
            spots[spotnr][4] = yl;
            spots[spotnr][5] = zl;
            spots[spotnr][6] = omegas[i];
            spots[spotnr][7] = etas[i];
            spots[spotnr][8] = theta;
            spots[spotnr][9] = ringnr;
            spotnr++;  // spot nr for this orientation             
            spotid++;  // overal spot id    
         }
      }
   }
   *nspots = spotnr;
}


 
////////////////////////////////////////////////////////////////////////////////
//
// returns the index of the minimum value in an array.
// 
// The array can be masked by array idxs
// 
void 
FindMinimum( RealType *array, int *idxs, int idxsSize, int *minIndex){ 
   int i;
   if (idxsSize == 0) { 
      *minIndex = -1;
      return;
   }
  
   *minIndex = idxs[0];
   for(i=0;i<idxsSize;++i)  { 
      if(array[idxs[i]] < array[*minIndex]) *minIndex = idxs[i];
   }
}

 
////////////////////////////////////////////////////////////////////////////////
//
// Compares a set of Theoretical spots with obs spots.
//
// returns the number of matches (within the given ranges)
// and for each match the difference between the theoretical and obs spot. 
// 
void    
CompareSpots(
   RealType **TheorSpots,
   int   nTheorSpots,
   RealType ObsSpots[][N_COL_OBSSPOTS],
//   int   nObsSpots,
   RealType RefRad,
   RealType MarginRad,
   RealType MarginRadial,
   RealType etamargins[],
   RealType omemargins[],
   int   *nMatch,
   RealType **GrainSpots)   
//   RealType **GrainSpots) 
   
{
   int nMatched = 0;      
   int nNonMatched = 0; 
   int sp;
   int RingNr;
   int iOme, iEta;
   int spotRow, spotRowBest;      
   int MatchFound ;
   RealType diffOme;
//   RealType diffRad;     
   RealType diffOmeBest;
   int iRing;
   int iSpot;
//   RealType diffIdRingDist;
   RealType etamargin, omemargin;
//   RealType diffEta;

   // for each spot in TheorSpots check if there is an equivalent one in ObsSpots
   for ( sp = 0 ; sp < nTheorSpots ; sp++ )  {
      RingNr = (int) TheorSpots[sp][9];
      iRing = RingNr-1;
      iEta = floor((180+TheorSpots[sp][12])/EtaBinSize);
      iOme = floor((180+TheorSpots[sp][6])/OmeBinSize);
      etamargin = etamargins[RingNr];
      omemargin = omemargins[(int) floor(fabs(TheorSpots[sp][12]))];  // omemargin depends on eta
      
#ifdef DEBUG            
   printf("iRing, iEta, iOme: %d %d %d\n", iRing, iEta, iOme );
#endif  
  
      // calc the diff. NOte: smallest in Omega difference is choosen
      MatchFound = 0;
      diffOmeBest = omemargin;
      for ( iSpot = 0 ; iSpot < ndata[iRing][iEta][iOme] ; iSpot++ ) {
         spotRow = data[iRing][iEta][iOme][iSpot];
         //diffIdRingDist = fabs(TheorSpots[sp][13] - ObsSpots[spotRow][8]);                  
         if ( fabs(TheorSpots[sp][13] - ObsSpots[spotRow][8]) < MarginRadial )  {
            //diffRad = fabs(RefRad - ObsSpots[spotRow][3]);                              
            if ( fabs(RefRad - ObsSpots[spotRow][3]) < MarginRad ) {            
               //diffEta = fabs(TheorSpots[sp][12] - ObsSpots[spotRow][6]);
               if ( fabs(TheorSpots[sp][12] - ObsSpots[spotRow][6]) < etamargin ) {
                  diffOme = fabs(TheorSpots[sp][6] - ObsSpots[spotRow][2]);                    
                  if ( diffOme < diffOmeBest ) {
                     diffOmeBest = diffOme;                  
                     spotRowBest = spotRow;
                     MatchFound = 1;
                  }
               }
            }
         }
      }
      
      if (MatchFound == 1) {
         // To be stored for each spot: 
         GrainSpots[nMatched][0] = nMatched;
         GrainSpots[nMatched][1] = 999.0;  // dummy
        
         GrainSpots[nMatched][2] = TheorSpots[sp][10];
         GrainSpots[nMatched][3] = ObsSpots[spotRowBest][0];
         GrainSpots[nMatched][4] = ObsSpots[spotRowBest][0] - TheorSpots[sp][10];       
         
         GrainSpots[nMatched][5] = TheorSpots[sp][11];
         GrainSpots[nMatched][6] = ObsSpots[spotRowBest][1];
         GrainSpots[nMatched][7] = ObsSpots[spotRowBest][1] - TheorSpots[sp][11];
        
         GrainSpots[nMatched][8] = TheorSpots[sp][6];
         GrainSpots[nMatched][9] = ObsSpots[spotRowBest][2];
         GrainSpots[nMatched][10]= ObsSpots[spotRowBest][2] - TheorSpots[sp][6];
        
         GrainSpots[nMatched][11] = RefRad;
         GrainSpots[nMatched][12] = ObsSpots[spotRowBest][3];
         GrainSpots[nMatched][13] = ObsSpots[spotRowBest][3] - RefRad;
         
         GrainSpots[nMatched][14] = ObsSpots[spotRowBest][4];
        
         nMatched++;    
      }
      else {  // the theoritcal spot is not found, store it at the end with negative id
         nNonMatched++;       
         int idx = nTheorSpots-nNonMatched;
         GrainSpots[idx][0] = -nNonMatched;
         GrainSpots[idx][1] = 999.0;  
         GrainSpots[idx][2] = TheorSpots[sp][10];
         GrainSpots[idx][3] = 0;         
         GrainSpots[idx][4] = 0;         
         GrainSpots[idx][5] = TheorSpots[sp][11];
         GrainSpots[idx][6] = 0;
         GrainSpots[idx][7] = 0;                  
         GrainSpots[idx][8] = TheorSpots[sp][6];
         GrainSpots[idx][9] = 0;                  
         GrainSpots[idx][10] = 0;
         GrainSpots[idx][11] = 0;                  
         GrainSpots[idx][12] = 0;         
         GrainSpots[idx][13] = 0;         
         GrainSpots[idx][14] = 0;         
      }
   }
  
   *nMatch = nMatched;  
}     
     
     


////////////////////////////////////////////////////////////////////////////////
// %%
// % Calculates the rotation matrix R from an axis/angle pair.
// %
// % NB: The rotation matrix R is such that: Crystal = R * Reference. 
// % Crystal is expressed in terms of the reference system!
// %
// % Input: 
// %    axis: a 3d vector (vector does not have to be normalized; this is done
// %          in the routine)
// %    angle: rotation angle [degrees]
// % 
// % Ouput: 3x3 rotation matrix.
// %
// % Date: 18-10-2010
// %
// % version 1.1 (11-1-2011): rearranged terms to make it a bit faster (~30%
// % faster)
// %
// 
// %%
void
AxisAngle2RotMatrix(
   RealType axis[3], 
   RealType angle,
   RealType R[3][3])  {
  
   // if axis is 0 0 0 then return just the Identy matrix
   if ( (axis[0] == 0) && (axis[1] == 0) && (axis[2] == 0) ) {
      R[0][0] = 1; 
      R[1][0] = 0;
      R[2][0] = 0;
      
      R[0][1] = 0;
      R[1][1] = 1;
      R[2][1] = 0;
      
      R[0][2] = 0;
      R[1][2] = 0;
      R[2][2] = 1;
      return;
   }
    
   // normalize vector first
   RealType lenInv = 1/sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
   RealType u = axis[0]*lenInv;
   RealType v = axis[1]*lenInv;    
   RealType w = axis[2]*lenInv;    
   RealType angleRad = deg2rad * angle;
   
   // source: http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q38
   // note: precalc of u*v and 1-rcos made it a bit slower (suprisingly)
   RealType rcos = cos(angleRad);
   RealType rsin = sin(angleRad);

   R[0][0] =      rcos + u*u*(1-rcos);
   R[1][0] =  w * rsin + v*u*(1-rcos);
   R[2][0] = -v * rsin + w*u*(1-rcos);
   
   R[0][1] = -w * rsin + u*v*(1-rcos);
   R[1][1] =      rcos + v*v*(1-rcos);
   R[2][1] =  u * rsin + w*v*(1-rcos);
   
   R[0][2] =  v * rsin + u*w*(1-rcos);
   R[1][2] = -u * rsin + v*w*(1-rcos);
   R[2][2] =      rcos + w*w*(1-rcos);
}
  	
  	
////////////////////////////////////////////////////////////////////////////////
// calculates the minimum angle to rotate around the normal of the give hkl plane, 
// too get an symetricaly equivalent orientation  
//  	if hkl are the same -> 120 degrees
//    if 2 are the same + 1 zero -> 180 
//    if 2 zeros -> 90
//     else 360
// NB: only works for cubic
//
RealType
RotationAngle(
  int hkl[3])
  
{
   int h, k, l;
   h = abs(hkl[0]);
   k = abs(hkl[1]);
   l = abs(hkl[2]);
   
   // count no of zeros
   int nzeros = 0;
   if ( h==0 ) nzeros++;
   if ( k==0 ) nzeros++;   
   if ( l==0 ) nzeros++;
   
   if (nzeros == 3) return 0;
   if (nzeros == 2) return 90;
   if (nzeros == 1) {  // check if the other 2 are equal
      if ( (h == 0) && ( k == l) ) return 180;
      if ( (k == 0) && ( h == l) ) return 180;
      if ( (l == 0) && ( h == k) ) return 180;
   }
   // check if hkl are all the same (and nonzero, but that is already cheched)
   if ( (h == k) && (k == l) ) return 120;
   
   // all other cases return 360
   return 360;
  
}
  
////////////////////////////////////////////////////////////////////////////////    
// % Calculate a list of unique orientations that have one plane in common (in this
// % case the diffraction plane).
// %
// % Input:
// %   hkl : the miller indices of the plane ie [ 1 1 0 ]
// %   hklnormal : the vector in 3d space that is the normal of the hkl plane
// %                (ie [ 1 3 2 ] )
// %   stepsize : stepsize of rotation around the plane normal. [degrees]
// %                           
// % Output:
// %   Orientations (rotation matrices) which all have the given hkl plane and
// %   hkl-normal (as input) in common. The 3x3 matrix is stored as 1x9 row
// %   (for each orientation 1 row).
// %
// % Date: oct 2010
// %
int
GenerateCandidateOrientationsF(
   int hkl[3], // miller indices of plane 
   RealType hklnormal[3],// direction of this plane in space (just the normal of the plane)
   RealType stepsize,
   RealType OrMat[][3][3],
   int * nOrient )  {
  
   RealType v[3];
   RealType MaxAngle = 0;   
  
   // calculate orientation via axis/angle pair
   crossProduct(v, hkl, hklnormal);      // axis is just the vector that is orthogonal to v1 and v2!
   RealType hkllen = sqrt(hkl[0]*hkl[0] + hkl[1]*hkl[1] + hkl[2]*hkl[2]);
   RealType hklnormallen = sqrt(hklnormal[0]*hklnormal[0] + hklnormal[1]*hklnormal[1] + hklnormal[2]*hklnormal[2]);
   RealType dotpr = dot(hkl, hklnormal); 
   RealType angled = rad2deg * acos(dotpr/(hkllen*hklnormallen));  // Angle is just the rotation in the plane defined by v1 v2!

   // calc the rotation matrix for this orientation (nb, one axis is assumed
   // now: vector v)
   // maybe faster way, using hkl and uvw
   RealType RotMat[3][3];
   RealType RotMat2[3][3];
   RealType RotMat3[3][3];
   AxisAngle2RotMatrix(v, angled, RotMat);
  
   // determine the rotation angle around the plane-normal (depends on the hkl plane)
   MaxAngle = RotationAngle(hkl);
   
#ifdef DEBUG
   printf("Rotation angle, maxangle: %lf\n", MaxAngle);
#endif    
  
   RealType nsteps = (MaxAngle/stepsize);
   int nstepsi = (int) nsteps;
   int or;
   int row, col;
   RealType angle2;
  
   // calculate all orientations 
   for ( or=0 ; or < nstepsi ; or++) {
      angle2 = or*stepsize;
      
      // now calculate the rotation matrix to rotate around the normal of the plane
      AxisAngle2RotMatrix(hklnormal, angle2, RotMat2);
      MatrixMultF33(RotMat2, RotMat, RotMat3);
      
      for (row = 0 ; row < 3 ; row++) {
         for (col = 0 ; col < 3 ; col++) {
             OrMat[or][row][col] = RotMat3[row][col];
         }
      }
   }
   *nOrient = nstepsi;
   
   return 0;      
}


////////////////////////////////////////////////////////////////////////////////
// %% 
// % Code to calculate displacement of spots needed if original position of
// % the center of mass of the grain is a,b,c when omega = 0. Vector in
// % direction of the diffracted beam is xi,yi,zi.
// %
// % All vectors (input and output) are in fable lab coordinates.
// %
// % Input: 
// %   a,b,c: position of grain [micron, mm, etc ]
// %   xi,yi,zi: vector of diffracted beam 
// %   omega: rotation of sample around the vertical axis (z), + is ccw [degrees]
// %
// % Output:
// %   Displ_y, displ_z: shift in y and z of the spot on 
// %        the detector.  [units same as input units]
// % 
// % Date: 22-10-2010
// %
void
displacement_spot_needed_COM(
   RealType a, 
   RealType b, 
   RealType c, 
   RealType xi, 
   RealType yi, 
   RealType zi, 
   RealType omega, 
   RealType *Displ_y, 
   RealType *Displ_z)
  
{
   RealType lenInv = 1/sqrt(xi*xi + yi*yi + zi*zi);
   xi = xi*lenInv;
   yi = yi*lenInv;
   zi = zi*lenInv;
  
   RealType OmegaRad = deg2rad * omega;
   RealType sinOme = sin(OmegaRad);
   RealType cosOme = cos(OmegaRad);
   RealType t = (a*cosOme - b*sinOme)/xi;
  
   *Displ_y = ((a*sinOme)+(b*cosOme)) -(t*yi);
   *Displ_z = c - t*zi;
}  


////////////////////////////////////////////////////////////////////////////////
// % Code to calculate g-vector from the spot position on the ring.
// % Important: The spot should be the idealized spot (from center of rotation) .
// %
// % INPUT:
// % - xi, yi, zi = vector in direction of the spot from the lab
// %                coordinate system. (diffracted ray)
// % - Omega [ degrees]: rotation angle of sample, when grain comes in diffraction.
////////////////////////////////////////////////////////////////////////////////
void
spot_to_gv(
   RealType xi,
   RealType yi,
   RealType zi,
   RealType Omega,
   RealType *g1,
   RealType *g2,
   RealType *g3)
{  
   // normalize
   RealType len = sqrt(xi*xi + yi*yi + zi*zi);
   
   if (len == 0) {
      *g1 = 0;
      *g2 = 0;
      *g3 = 0;
      printf("len o!\n");            
      return;
   }
   
   RealType xn = xi/len;
   RealType yn = yi/len;
   RealType zn = zi/len;

   RealType g1r = (-1 + xn);  // r means rotated (gvec at rotation omega)
   RealType g2r = yn;
  
   //% rotate back to omega = 0
   RealType CosOme = cos(-Omega*deg2rad);
   RealType SinOme = sin(-Omega*deg2rad);   

   *g1 = g1r * CosOme - g2r * SinOme;
   *g2 = g1r * SinOme + g2r * CosOme;
   *g3 = zn;   
   // g3 (z) does not change during rotation   
}


////////////////////////////////////////////////////////////////////////////////
// % Code to calculate g-vector from the spot position on the ring.
// % This version takes into account the pos of the grain in the sample.
// %
// % INPUT:
// % - xi, yi, zi = vector in direction of the spot from the lab
// %                coordinate system.
// % - Omega [ degrees]
// % - cx, cy, cz = pos of the crystal (at omega = 0)
// %
// % OUTPUT:
// % - g vector (at omega = 0 !)
// %
////////////////////////////////////////////////////////////////////////////////
void
spot_to_gv_pos(
   RealType xi,
   RealType yi,
   RealType zi,
   RealType Omega,
   RealType cx,
   RealType cy,
   RealType cz,
   RealType *g1,
   RealType *g2,
   RealType *g3)
   
{  
   RealType v[3], vr[3];
   
   // first correct the vector xi, yi, zi for the pos of grain in the sample.
   // subtract the grain pos (use pos when it is rotated! not the normal pos (ome = 0)
   v[0] = cx;
   v[1] = cy;
   v[2] = cz;
   RotateAroundZ(v, Omega, vr); // vr is the pos of the grain after rotation
   xi = xi - vr[0];
   yi = yi - vr[1];
   zi = zi - vr[2];
   
   spot_to_gv( xi, yi, zi, Omega, g1, g2, g3);
}


////////////////////////////////////////////////////////////////////////////////
void
FriedelEtaCalculation(
   RealType ys,
   RealType zs,
   RealType ttheta,
   RealType eta,
   RealType Ring_rad,
   RealType Rsample,
   RealType Hbeam, 
   RealType *EtaMinFr,
   RealType *EtaMaxFr)
   
{
   RealType quadr_coeff2 = 0;
   RealType eta_Hbeam, quadr_coeff, coeff_y0 = 0, coeff_z0 = 0, y0_max_z0, y0_min_z0, y0_max = 0, y0_min = 0, z0_min = 0, z0_max = 0;

   // % Calculate some parameters (These can also be calculated once for each hkl
   // % set and then called during the program).
   
   //RealType Ring_rad = Lsd*tan(ttheta*deg2rad);
   
   if (eta > 90)
     eta_Hbeam = 180 - eta;
   else if (eta < -90)
     eta_Hbeam = 180 - fabs(eta);
   else
     eta_Hbeam = 90 - fabs(eta);
   
   Hbeam = Hbeam + 2*(Rsample*tan(ttheta*deg2rad))*(sin(eta_Hbeam*deg2rad));
   
   RealType eta_pole = 1 + rad2deg*acos(1-(Hbeam/Ring_rad));
   RealType eta_equator = 1 + rad2deg*acos(1-(Rsample/Ring_rad));
   
   // Find out which quadrant is the spot in
   if ((eta >= eta_pole) && (eta <= (90-eta_equator)) ) { // % 1st quadrant
      quadr_coeff = 1;
      coeff_y0 = -1;                                      
      coeff_z0 = 1;
   }
   else if ( (eta >=(90+eta_equator)) && (eta <= (180-eta_pole)) ) {//% 4th quadrant
      quadr_coeff = 2;
      coeff_y0 = -1;
      coeff_z0 = -1;
   }
   else if ( (eta >= (-90+eta_equator) ) && (eta <= -eta_pole) )   { // % 2nd quadrant
      quadr_coeff = 2;
      coeff_y0 = 1;
      coeff_z0 = 1;
   }  
   else if ( (eta >= (-180+eta_pole) ) && (eta <= (-90-eta_equator)) )  { // % 3rd quadrant
      quadr_coeff = 1;
      coeff_y0 = 1;
      coeff_z0 = -1;
   }
   else
     quadr_coeff = 0;
   
   
   //% Calculate y0 max and min due to Rsample.
   RealType y0_max_Rsample = ys + Rsample;
   RealType y0_min_Rsample = ys - Rsample;
   
   // Calculate z0 max and min due to Hbeam
   RealType z0_max_Hbeam = zs + 0.5 * Hbeam;
   RealType z0_min_Hbeam = zs - 0.5 * Hbeam;
   
   // Calculate y0 max and min due to z0
   if (quadr_coeff == 1) {
      y0_max_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_max_Hbeam * z0_max_Hbeam));
      y0_min_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_min_Hbeam * z0_min_Hbeam));
   }
   else if (quadr_coeff == 2) {
      y0_max_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_min_Hbeam * z0_min_Hbeam));
      y0_min_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_max_Hbeam * z0_max_Hbeam));
   }
   
   
   // Select whether to choose the limit due to Rsample or due to z0
   if (quadr_coeff > 0)  {
      y0_max = min(y0_max_Rsample, y0_max_z0);
      y0_min = max(y0_min_Rsample, y0_min_z0);
   }
   else {
      if ((eta > -eta_pole) && (eta < eta_pole ))  {
          y0_max = y0_max_Rsample;
          y0_min = y0_min_Rsample;
          coeff_z0 = 1;
      }
      else if (eta < (-180+eta_pole))  {
          y0_max = y0_max_Rsample;
          y0_min = y0_min_Rsample;
          coeff_z0 = -1;
      }
      else if (eta > (180-eta_pole))  {
          y0_max = y0_max_Rsample;
          y0_min = y0_min_Rsample;
          coeff_z0 = -1;
      }
      else if (( eta > (90-eta_equator)) && (eta < (90+eta_equator)) ) {
          quadr_coeff2 = 1;
          z0_max = z0_max_Hbeam;
          z0_min = z0_min_Hbeam;
          coeff_y0 = -1;
      }
      else if ((eta > (-90-eta_equator)) && (eta < (-90+eta_equator)) ) {
          quadr_coeff2 = 1;
          z0_max = z0_max_Hbeam;
          z0_min = z0_min_Hbeam;
          coeff_y0 = 1;
      }
   }
   
   // Calculate y0_min (max) or z0_min (max).
   if ( quadr_coeff2 == 0 ) {
       z0_min = coeff_z0 * sqrt((Ring_rad * Ring_rad)-(y0_min * y0_min));
       z0_max = coeff_z0 * sqrt((Ring_rad * Ring_rad)-(y0_max * y0_max));
   }
   else {
       y0_min = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_min * z0_min));
       y0_max = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_max * z0_max));
   }
   
   RealType dYMin = ys - y0_min;
   RealType dYMax = ys - y0_max;
   RealType dZMin = zs - z0_min;
   RealType dZMax = zs - z0_max;
   
   // Calculate the Friedel pair locations on ideal ring
   RealType YMinFrIdeal =  y0_min;
   RealType YMaxFrIdeal =  y0_max;
   RealType ZMinFrIdeal = -z0_min;
   RealType ZMaxFrIdeal = -z0_max;
   
   // Calculate the Extremities on the displaced location
   RealType YMinFr = YMinFrIdeal - dYMin;
   RealType YMaxFr = YMaxFrIdeal - dYMax;
   RealType ZMinFr = ZMinFrIdeal + dZMin;
   RealType ZMaxFr = ZMaxFrIdeal + dZMax;
   
   // Calculate etas for friedel circle
   RealType Eta1, Eta2;
   CalcEtaAngle((YMinFr + ys),(ZMinFr - zs), &Eta1);
   CalcEtaAngle((YMaxFr + ys),(ZMaxFr - zs), &Eta2);
   
   *EtaMinFr = min(Eta1,Eta2);
   *EtaMaxFr = max(Eta1,Eta2);
}   


RealType
sign(RealType anumber) 
{
  if (anumber < 0) return -1.0;
  else return 1.0;

}

////////////////////////////////////////////////////////////////////////////////
// %%
// % Code to step in orientation by changing the location of the ideal spot
// % along the ring. 
// %
// % Input:
// %  - ys, zs: lab coordinates of real spot [micron]
// %  - 2theta: diffraction angle [degrees]
// %  - eta: azimuth angle + = ccw [degrees]
// %  - Ring_rad: radius of the diffraction ring [micron]
// %  - Rsample : radius of sample [micron]
// %  - Hbeam: Height of the beam [micron]
// %  - step_size: Step size in the sample [microns]
// %
// % Output: 
// %  - y0_vector, z0_vector: the coordinates of ideal spots for a grain in the 
// %    center of the sample (lab coordinates). [micron]
// %
// 

void
GenerateIdealSpots(
   RealType ys, 
   RealType zs, 
   RealType ttheta, 
   RealType eta, 
   RealType Ring_rad, 
   RealType Rsample, 
   RealType Hbeam, 
   RealType step_size,
   RealType y0_vector[], 
   RealType z0_vector[], 
   int * NoOfSteps) 
  
{
   int quadr_coeff2 = 0;
   RealType eta_Hbeam, quadr_coeff, coeff_y0 = 0, coeff_z0 = 0, y0_max_z0, y0_min_z0, y0_max = 0, y0_min = 0, z0_min = 0, z0_max = 0;
   RealType y01, z01, y02, z02, y_diff, z_diff, length;
   int nsteps;  
   
   // % Calculate some parameters (These can also be calculated once for each hkl
   // % set and then called during the program).
   
   //RealType Ring_rad = Lsd*tan(ttheta*deg2rad);
   
   if (eta > 90)
     eta_Hbeam = 180 - eta;
   else if (eta < -90)
     eta_Hbeam = 180 - fabs(eta);
   else
     eta_Hbeam = 90 - fabs(eta);
   
   Hbeam = Hbeam + 2*(Rsample*tan(ttheta*deg2rad))*(sin(eta_Hbeam*deg2rad));
   
   RealType eta_pole = 1 + rad2deg*acos(1-(Hbeam/Ring_rad));
   RealType eta_equator = 1 + rad2deg*acos(1-(Rsample/Ring_rad));
   
   // Find out which quadrant is the spot in
   if ((eta >= eta_pole) && (eta <= (90-eta_equator)) ) { // % 1st quadrant
      quadr_coeff = 1;
      coeff_y0 = -1;                                      
      coeff_z0 = 1;
   }
   else if ( (eta >=(90+eta_equator)) && (eta <= (180-eta_pole)) ) {//% 4th quadrant
      quadr_coeff = 2;
      coeff_y0 = -1;
      coeff_z0 = -1;
   }
   else if ( (eta >= (-90+eta_equator) ) && (eta <= -eta_pole) )   { // % 2nd quadrant
      quadr_coeff = 2;
      coeff_y0 = 1;
      coeff_z0 = 1;
   }  
   else if ( (eta >= (-180+eta_pole) ) && (eta <= (-90-eta_equator)) )  { // % 3rd quadrant
      quadr_coeff = 1;
      coeff_y0 = 1;
      coeff_z0 = -1;
   }
   else
     quadr_coeff = 0;
   
   
   //% Calculate y0 max and min due to Rsample.
   RealType y0_max_Rsample = ys + Rsample;
   RealType y0_min_Rsample = ys - Rsample;
   
   // Calculate z0 max and min due to Hbeam
   RealType z0_max_Hbeam = zs + 0.5 * Hbeam;
   RealType z0_min_Hbeam = zs - 0.5 * Hbeam;
   
   // Calculate y0 max and min due to z0
   if (quadr_coeff == 1) {
      y0_max_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_max_Hbeam * z0_max_Hbeam));
      y0_min_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_min_Hbeam * z0_min_Hbeam));
   }
   else if (quadr_coeff == 2) {
      y0_max_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_min_Hbeam * z0_min_Hbeam));
      y0_min_z0 = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_max_Hbeam * z0_max_Hbeam));
   }
   
   
   // Select whether to choose the limit due to Rsample or due to z0
   if (quadr_coeff > 0)  {
      y0_max = min(y0_max_Rsample, y0_max_z0);
      y0_min = max(y0_min_Rsample, y0_min_z0);
   }
   else {
      if ((eta > -eta_pole) && (eta < eta_pole ))  {
          y0_max = y0_max_Rsample;
          y0_min = y0_min_Rsample;
          coeff_z0 = 1;
      }
      else if (eta < (-180+eta_pole))  {
          y0_max = y0_max_Rsample;
          y0_min = y0_min_Rsample;
          coeff_z0 = -1;
      }
      else if (eta > (180-eta_pole))  {
          y0_max = y0_max_Rsample;
          y0_min = y0_min_Rsample;
          coeff_z0 = -1;
      }
      else if (( eta > (90-eta_equator)) && (eta < (90+eta_equator)) ) {
          quadr_coeff2 = 1;
          z0_max = z0_max_Hbeam;
          z0_min = z0_min_Hbeam;
          coeff_y0 = -1;
      }
      else if ((eta > (-90-eta_equator)) && (eta < (-90+eta_equator)) ) {
          quadr_coeff2 = 1;
          z0_max = z0_max_Hbeam;
          z0_min = z0_min_Hbeam;
          coeff_y0 = 1;
      }
   }
   
   //% calculate nsteps
   if (quadr_coeff2 == 0 ) {
      y01 = y0_min;
      z01 = coeff_z0 * sqrt((Ring_rad * Ring_rad )-(y01 * y01));
      y02 = y0_max;
      z02 = coeff_z0 * sqrt((Ring_rad * Ring_rad )-(y02 * y02));
      y_diff = y01 - y02;
      z_diff = z01 - z02;
      length = sqrt(y_diff * y_diff + z_diff * z_diff);
      nsteps = ceil(length/step_size);
   }
   else {
      z01 = z0_min;
      y01 = coeff_y0 * sqrt((Ring_rad * Ring_rad )-((z01 * z01)));
      z02 = z0_max;
      y02 = coeff_y0 * sqrt((Ring_rad * Ring_rad )-((z02 * z02)));
      y_diff = y01 - y02;
      z_diff = z01 - z02;
      length = sqrt(y_diff * y_diff + z_diff * z_diff);
      //nsteps = (int)length/step_size;
      nsteps = ceil(length/step_size);      
   }
   
   // make nsteps odd, to make sure we have a middle point
   if ((nsteps % 2) == 0 ) {
     nsteps = nsteps +1;
   }
   
   // special case: if no of steps is only 1 take the middle point
   if ( nsteps == 1 ) {
      if (quadr_coeff2 == 0) {
         y0_vector[0] = (y0_max+y0_min)/2;
         z0_vector[0] = coeff_z0 * sqrt((Ring_rad * Ring_rad)-(y0_vector[0] * y0_vector[0]));
      }
      else {
         z0_vector[0] = (z0_max+z0_min)/2;
         y0_vector[0] = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_vector[0] * z0_vector[0]));
      }
   }
   else {
      int i;
      RealType stepsizeY = (y0_max-y0_min)/(nsteps-1);
      RealType stepsizeZ = (z0_max-z0_min)/(nsteps-1);
      
      
      if (quadr_coeff2 == 0) {
         for (i=0 ; i < nsteps ; i++) {
            y0_vector[i] = y0_min + i*stepsizeY;
            z0_vector[i] = coeff_z0 * sqrt((Ring_rad * Ring_rad)-(y0_vector[i] * y0_vector[i]));
         }
      }
      else {
         for (i=0 ; i < nsteps ; i++) {  
            z0_vector[i] = z0_min + i*stepsizeZ;
            y0_vector[i] = coeff_y0 * sqrt((Ring_rad * Ring_rad)-(z0_vector[i] * z0_vector[i]));
           
         }
      }
   }
   *NoOfSteps = nsteps;        
}



////////////////////////////////////////////////////////////////////////////////   
// % Code to calculate n_max and n_min for use in the code
// % spot_to_unrotated_coordinates (should be run after this code). It needs
// % to be run for each combination of ys and y0 as obtained from the 
// % stepping_orientation_spot code.
// %
// % INPUT
// % - xi,yi: Unit vector in direction of the diffracted beam. [microns]
// % - ys, y0: y-position of the spot in Lab coordinates for the real spot and
// %           ideal spot along the ring. [microns]
// % - Rsample: radius of the sample [microns]
// % - step_size: step size in the sample [microns]
// %
// % OUTPUT
// % - n_max, n_min: Inputs in the spot_to_unrotated_coordinates code.
//
//

void
calc_n_max_min(
   RealType xi, 
   RealType yi, 
   RealType ys, 
   RealType y0, 
   RealType R_sample, 
   int   step_size, 
   int * n_max, 
   int * n_min)
{
   RealType dy = ys-y0;
   RealType a = xi*xi + yi*yi;
   RealType b = 2*yi*dy;
   RealType c = dy*dy - R_sample*R_sample;
   RealType D = b*b - 4*a*c;
//   if (D < 0) printf("Error calculating number of steps in sample\n");
   RealType P = sqrt(D);
   RealType lambda_max = (-b+P)/(2*a) + 20; // 20 microns is the inaccuracy in aligning the sample.

   //% Calculate the values of n_max and n_min
   *n_max = (int)((lambda_max*xi)/(step_size));
   *n_min = - *n_max;
}




//////////////////////////////////////////////////////////////////////////////// 
// %%
// % Code to find the un-rotated coordinates (a,b,c) of the center of mass of
// % a grain from a spot. Assumes an orientation along the ring (which can be 
// % stepped through using the code stepping_orientation_spot), and a plane in 
// % the sample along the beam direction (this can also be stepped through using
// % lambda equation given below). 
// %
// % Input:
// %   - xi,yi,zi: Unit vector along the diffracted beam.  
// %   - ys,zs: Spot coordinates on detector (normalized w.r.t. beam center).
// %            [micron]
// %   - y0,z0: coordinates according to orientation at the ring
// %            (normalized w.r.t. beam center). [units same as ys and zs]
// %   - step_size_in_x: the step in x [micron]
// %   - n : goes from -RadiusOfSample/stepsize to +RadiusOfSample/stepsize.
// %        0 = center of sample.
// %   - omega: rotation of sample around vertical axis + = ccw [degrees]
// %
// % Output:
// %   a, b, c: unrotated coordinates of the (center of mass) of the grain. 
// %           [units same as ys and zs] 
// %
// % Date: 22-10-2010
// %

void
spot_to_unrotated_coordinates(
   RealType xi, 
   RealType yi, 
   RealType zi, 
   RealType ys, 
   RealType zs, 
   RealType y0,  
   RealType z0, 
   RealType step_size_in_x, 
   int   n, 
   RealType omega,
   RealType *a,
   RealType *b,
   RealType *c)
{
   RealType lambda = (step_size_in_x)*(n/xi);
   //   % Where n can be from -100 to 100 in case of a 1 mm thick sample and 5
   //   % microns step size.

   //% Rotated coordinates
   RealType x1 = lambda*xi;
   RealType y1 = ys - y0 + lambda*yi;
   RealType z1 = zs - z0 + lambda*zi;
  
   //% Un-rotated coordinates
   RealType cosOme = cos(omega*deg2rad);
   RealType sinOme = sin(omega*deg2rad);
   *a = (x1*cosOme) + (y1*sinOme);
   *b = (y1*cosOme) - (x1*sinOme);
   *c = z1;
}




////////////////////////////////////////////////////////////////////////////////
// tries to find the friedel pair. And returns coordinates (y,z) of ideal spots
// on the ring (which defines the diffracted rays, which gives the planenormals). 
//
//
void
GenerateIdealSpotsFriedel(
   RealType ys, 
   RealType zs, 
   RealType ttheta, 
   RealType eta, 
   RealType omega,
   int ringno,
   RealType Ring_rad, 
   RealType Rsample, 
   RealType Hbeam,
   RealType OmeTol, // tolerance for omega difference of the friedelpair
   RealType RadiusTol,    // tolerance for difference of the friedelpair in radial direction 
   RealType y0_vector[],   // returns: y and z coordinates of the spots found
   RealType z0_vector[], 
   int * NoOfSteps)   // no of spots found
    
{
   RealType EtaF;
   RealType OmeF;
   RealType EtaMinF, EtaMaxF, etaIdealF;
   RealType IdealYPos, IdealZPos;
   *NoOfSteps = 0;
   
   if (omega < 0 )  OmeF = omega + 180;
   else             OmeF = omega - 180;
    
   if ( eta < 0 )  EtaF = -180 - eta;
   else            EtaF = 180 - eta;
   
   // now find the candidate friedel spots
   int r;
   int rno_obs;
   RealType ome_obs, eta_obs;
   //printf("omega, OmeFr: %lf %lf\n", omega, OmeF);
   for (r=0 ; r < n_spots ; r++) {
       // filter on ringno
       rno_obs = round(ObsSpotsLab[r][5]);
       ome_obs = ObsSpotsLab[r][2];
       eta_obs = ObsSpotsLab[r][6];
       
       // skip if outside margins
       if ( rno_obs != ringno ) continue;
       if ( fabs(ome_obs - OmeF) > OmeTol) continue;
       
       // radial margins
       RealType yf = ObsSpotsLab[r][0];
       RealType zf = ObsSpotsLab[r][1];
       RealType EtaTransf;  // NB this eta is defined differently: the orignal spot is the origin now.
       CalcEtaAngle(yf + ys, zf - zs, &EtaTransf); 
       RealType radius = sqrt((yf + ys)*(yf + ys) + (zf - zs)*(zf - zs));
       if ( fabs(radius - 2*Ring_rad) > RadiusTol)  continue;
       
       // calculate eta boundaries
       FriedelEtaCalculation(ys, zs, ttheta, eta, Ring_rad, Rsample, Hbeam, &EtaMinF, &EtaMaxF);
              
       // in eta segment?
       if (( EtaTransf < EtaMinF) || (EtaTransf > EtaMaxF) ) continue;
       
       // found a candidate friedel pair!
       // calculate ideal spot (on the ring-> gives orientation)
       RealType ZPositionAccZ = zs - (( zf + zs)/2);
       RealType YPositionAccY = ys - ((-yf + ys)/2);

       CalcEtaAngle(YPositionAccY, ZPositionAccZ, &etaIdealF);
       CalcSpotPosition(Ring_rad, etaIdealF, &IdealYPos, &IdealZPos); 

       // save spot positions                               
       y0_vector[*NoOfSteps] = IdealYPos;
       z0_vector[*NoOfSteps] = IdealZPos;       
       (*NoOfSteps)++;
   }
}


int
AddUnique(int *arr, int *n, int val)
{
   int i;
   
   for (i=0 ; i < *n ; ++i) {
      if (arr[i] == val) {
          return 0;  // dont add
      }
   }
   
   // not in array, so add
   arr[*n] = val;
   (*n)++;
   return 1;
}



void
MakeUnitLength(
   RealType x,
   RealType y,
   RealType z, 
   RealType *xu, 
   RealType *yu, 
   RealType *zu )
{
   RealType len = CalcLength(x, y, z);
   
   if (len == 0) {
      *xu = 0;
      *yu = 0;
      *zu = 0;
      return; 
   }
   
   *xu = x/len;
   *yu = y/len;
   *zu = z/len;
}


////////////////////////////////////////////////////////////////////////////////
// tries to find the mixed friedel pair of a plane hkl (-h-k-l) And returns 
// coordinates (y,z) of ideal spots on the ring (which defines the 
// diffracted rays, which gives the planenormals).
// 
//%
//%% IMPORTANT: Please note that this works for spots with asind(abs(sind(eta))) > 10 degrees only.
//%
//%  Version: 1
//
void
GenerateIdealSpotsFriedelMixed(
   RealType ys,          // y coord of starting spot
   RealType zs,          // z coord ,,      
   RealType Ttheta,      // 2-theta ,,
   RealType Eta,         // eta ,,
   RealType Omega,       // omega ,,
   int RingNr,           // ring number ,,
   RealType Ring_rad,    // ring radius of ideal ring [micron]
   RealType Lsd,         // sample detector distance [micron]
   RealType Rsample,     // radius of sample [micron]
   RealType Hbeam,       // height of beam [micron]
   RealType StepSizePos, // stepsize of plane normals (and along the beam) [micron]   
   RealType OmeTol,      // tolerance for omega difference of the friedelpair
   RealType RadialTol,   // tolerance for difference of the friedelpair in radial direction
   RealType EtaTol,      // tolerance for eta difference of the fp  
   RealType spots_y[],   // returns: y and z coordinates of the ideal spots that gave a fp 
   RealType spots_z[], 
   int * NoOfSteps)      // no of ideal spots found
   
{
   const int MinEtaReject = 10; // [degrees] dont try spots near the pole.

   RealType omegasFP[4];
   RealType etasFP[4];
   int nsol;
   int nFPCandidates;
   RealType theta = Ttheta/2;
   RealType SinMinEtaReject = sin(MinEtaReject * deg2rad);
   RealType y0_vector[2000]; // yz of ideal spots
   RealType z0_vector[2000];
   RealType G1, G2, G3;
   int SpOnRing, NoOfSpots;
   int FPCandidatesUnique[2000]; // just stores the id's of the planenormals 
   RealType FPCandidates[2000][3];
   RealType xi, yi, zi;
   RealType y0, z0;
   RealType YFP1, ZFP1;
   int nMax, nMin;

   RealType EtaTolDeg;
   
   
   EtaTolDeg = rad2deg* atan(EtaTol / Ring_rad) ;// convert from micron to degrees
    
   // some init
   *NoOfSteps = 0;
   nFPCandidates = 0;

   // check border situation
   if (fabs(sin(Eta * deg2rad)) < SinMinEtaReject) {
      printf("The spot is too close to the poles. This technique to find mixed friedel pair would not work satisfactorily. So don't use mixed friedel pair.\n");
      return;
   }

   // Try points on the 'surface':
   // - calc ideal spots on ring: these are the plane normals 
   // - for each plane normal walk along the ray
   // - calc ideal spot pos of mixed fp
   // - displace spot due to pos in sample
   // - and check if this spot is in the dataset: if yes potential fp spot.
   // in the end returns for each hit, the best planenormals (=smallest difference in pos) 

   // first generate ideal spots
   GenerateIdealSpots(ys, zs, Ttheta, Eta, Ring_rad, Rsample, Hbeam, StepSizePos, y0_vector, z0_vector, &NoOfSpots);

   // check each planenormal    
   for (SpOnRing = 0 ; SpOnRing < NoOfSpots ; ++SpOnRing ) {
//       printf("spotno of spot: %d %d", SpOnRing, NoOfSpots);
       y0 = y0_vector[SpOnRing];
       z0 = z0_vector[SpOnRing];

       // unit vector in the direction of the ideal spot.  
       MakeUnitLength(Lsd, y0, z0, &xi, &yi, &zi);

       // calc gv and when the back of it comes into diffraction (omega)
       spot_to_gv(Lsd, y0, z0, Omega, &G1, &G2, &G3);
       CalcOmega(-G1, -G2, -G3, theta, omegasFP, etasFP, &nsol);     // take the back of plane!  
        
       // if no solutions go to next planenormal (if 1 solution: then there is also no mixed fp)
       if (nsol <= 1) {
          printf("no omega solutions. skipping plane.\n");
          continue;
       }
        
       // take the FP Omega that is closest to the original omega (this is the mixed fp)
       RealType OmegaFP, EtaFP, diff0, diff1;
       diff0 = fabs(omegasFP[0] - Omega);
       if (diff0 > 180)  { diff0 = 360 - diff0;  }  // in case the omegas are close to 180 and -180 the difference is big

       diff1 = fabs(omegasFP[1] - Omega);
       if (diff1 > 180)  { diff1 = 360 - diff1;  }  // in case the omegas are close to 180 and -180 the difference is big 

       // use the smallest
       if (  diff0 < diff1)  {
          OmegaFP = omegasFP[0];
          EtaFP   = etasFP[0];  
       }     
       else {
          OmegaFP = omegasFP[1];
          EtaFP   = etasFP[1];        
       }
        
       // calc the (ideal) spot pos of this mixed fp.
       CalcSpotPosition(Ring_rad, EtaFP, &YFP1, &ZFP1);
        
       // now step along the ray
       calc_n_max_min(xi, yi, ys, y0, Rsample, StepSizePos, &nMax, &nMin);
//       printf("xi, yi, ys, y0, Rsample, StepSizePos, RadialTol, &nMax, &nMin: %f %f %f %f %f %f %f %d %d\n", xi, yi, ys, y0, Rsample, StepSizePos, RadialTol, nMax, nMin);       
       RealType a, b, c, YFP, ZFP, RadialPosFP, EtaFPCorr;
       int n;
       for (n = nMin ; n <= nMax ; ++n ){
           spot_to_unrotated_coordinates(xi, yi, zi, ys, zs, y0, z0, StepSizePos, n, Omega, &a, &b, &c); 
               
            // if outside sample go to next n
            if (fabs(c) > Hbeam/2) {   continue;    }
            
            // now displace spot due to displacement of grain in sample
            RealType Dy,Dz;
            displacement_spot_needed_COM(a, b, c, Lsd, YFP1, ZFP1, OmegaFP, &Dy, &Dz); 
            YFP = YFP1 + Dy;
            ZFP = ZFP1 + Dz;
            
            RadialPosFP = sqrt(YFP*YFP + ZFP*ZFP) - Ring_rad;
            CalcEtaAngle(YFP,ZFP, &EtaFPCorr);
            
            // now check if this spot is in the dataset (using the bins only, a bit crude but fast way )
            int *spotRows;
            int nspotRows, iSpot, spotRow;            
            GetBin(RingNr, EtaFPCorr, OmegaFP, &spotRows, &nspotRows);
           
            RealType diffPos2, dy, dz;
            
            // check for each spot in the bin, if it close enough to the calculated spot
            for ( iSpot = 0 ; iSpot < nspotRows ; iSpot++ ) {
                spotRow = spotRows[iSpot];
                if ( (fabs(RadialPosFP - ObsSpotsLab[spotRow][8]) < RadialTol ) &&
                     (fabs(OmegaFP - ObsSpotsLab[spotRow][2]) < OmeTol ) &&                
                     (fabs(EtaFPCorr - ObsSpotsLab[spotRow][6]) < EtaTolDeg )  )
                {
                    // found a candidate fp: add (if not in the list already);
                    // if already in the list: keep the one with smallest diff
                    dy = (YFP-ObsSpotsLab[spotRow][0]);
                    dz = (ZFP-ObsSpotsLab[spotRow][1]);
                    diffPos2 = dy*dy + dz*dz;
                                               
                    int i;
                    int idx = nFPCandidates; // defaults to add 
                    for (i=0 ; i < nFPCandidates ; ++i) {
                        // check spot id
                        if (FPCandidates[i][0] == ObsSpotsLab[spotRow][4] )   {
                           // already in list: check if the diff is smaller
//                           printf("spotRow diffPos2 FPCandidates: %d %f %f\n", spotRow, diffPos2 ,FPCandidates[i][2]);
                           if (diffPos2 < FPCandidates[i][2] )  {
                              idx = i;  // item i will be replaced now
                           } else {
                              idx = -1; // indicating nothing should be done
                           }
//                            printf("diffPos2 stored diffr rtol diffw wtol diffe etol: %f %f %f %f %f %f %f %f\n", diffPos2 , FPCandidates[i][2],
//                      fabs(RadialPosFP - ObsSpotsLab[spotRow][8]),RadialTol ,
//                      fabs(OmegaFP - ObsSpotsLab[spotRow][2]) , OmeTol ,                
//                      fabs(EtaFPCorr - ObsSpotsLab[spotRow][6]) , EtaTolDeg )  ;

                           
                           break;
                           
                        }
                    }

                    // add or replace the new FP candidate  
                    // table FPCandidates: spotid_of_FP sponring diffpos                    
                    if (idx >= 0) {
                       FPCandidates[idx][0] = ObsSpotsLab[spotRow][4]; 
                       FPCandidates[idx][1] = SpOnRing;
                       FPCandidates[idx][2] = diffPos2;                 
                       if (idx == nFPCandidates) nFPCandidates++;
                    }                                            
                }
            }
        }
    }
    
    int i;
    // finally remove any double planenormals (a pn can be found multiple times for different spots)
    int nFPCandidatesUniq = 0;
    for (i=0 ; i < nFPCandidates ; ++i) {
       AddUnique(FPCandidatesUnique, &nFPCandidatesUniq, FPCandidates[i][1]);
    }
    
    // return the spot pos      
    int iFP;    
    for (iFP = 0 ; iFP < nFPCandidatesUniq ; ++iFP) {
       spots_y[iFP] = y0_vector[FPCandidatesUnique[iFP]];
       spots_z[iFP] = z0_vector[FPCandidatesUnique[iFP]];
//       printf("%d %d %f %f \n", iFP, FPCandidatesUnique[iFP], spots_y[iFP], spots_z[iFP]);
    }
    *NoOfSteps = nFPCandidatesUniq;      
}      



////////////////////////////////////////////////////////////////////////////////
// % Calculates diffraction angle theta 
// %
// % Input:
// %  - hkl plane (at least one should be nonzero)
// %  - LatticeParameter [A] 
// %  - Wavelength [A]
// %
// % Output: 
// %  Theta [degrees]
// %
// % Date: Nov 2010
// %
void
CalcTheta(
   int h, 
   int k, 
   int l, 
   RealType LatticeParameter, 
   RealType Wavelength, 
   RealType *theta) 
{
   RealType dspacing; 

   RealType h2k2l2 = h*h + k*k + l*l;
   dspacing = sqrt(LatticeParameter * LatticeParameter/h2k2l2);
   *theta = rad2deg*asin(Wavelength/(2*dspacing));
}



  


////////////////////////////////////////////////////////////////////////////////
// Calculates 2-theta [deg] from spotpos (y,z) [micron] and distance [micron]
//
void
CalcTtheta(
   RealType y, 
   RealType z,
   RealType distance,
   RealType *ttheta)
{
   RealType rad = sqrt((y*y) + (z*z));
   *ttheta = rad2deg * atan(rad/distance);
}

    
////////////////////////////////////////////////////////////////////////////////
// Parameters in the parameter file
//    
struct TParams {
   int RingNumbers[MAX_N_RINGS];   // the ring numbers to use for indexing (1, 2, 4, etc) 
   int CellStruct;                 // 1=bcc 2=fcc
   RealType LatticeConstant;          // [Angstrom]
   RealType Wavelength;               // Wavelength of incoming beam [Angstrom] 
   RealType Distance;                 // Distance between sample and detector [micron]
   RealType Rsample;                  // Radius of the sample [micron] 
   RealType Hbeam;                    // Height of the beam [micron]    
//   char HKLsFileName[1024];       
   RealType StepsizePos;              // step size in position [micron]  
   RealType StepsizeOrient;           // step size in orientation (rotation around the plane normal) [degrees]
   int NrOfRings;                  // No of rings to use (not explicit input by user, but set via RingNumbers[])
   RealType RingRadii[MAX_N_RINGS];   // Radii of the rings [micron]. this is a used internally: ringrad of ring 1 is at index 1 etc. 
   RealType RingRadiiUser[MAX_N_RINGS];   // Radii of the rings [micron]. stores only radii of the used rings!! Used for user input.    
   RealType MarginOme;                // Margin in Omega [degrees], when assigning theoretical spots to experimental spots. (|omeT-omeO| < MarginOme)
   RealType MarginEta;                // Margin in eta [degrees], ,,
   RealType MarginRad;                // Margin in radius [micron], ,, 
   RealType MarginRadial;             // Margin in radial direction (ortogonal to the ring) [micron], ,,
   RealType EtaBinSize;               // Size of bin for eta [degrees]
   RealType OmeBinSize;               // Size of bin for omega [degrees]
   RealType ExcludePoleAngle;         // Spots can be excluded at the poles: the range is |Eta| < ExcludePoleAngle and 180-|Eta| < ExcludePoleAngle [degrees]
   RealType MinMatchesToAcceptFrac;   // Minimum fraction (matched_spots/exp_spots) to accept an orientation+position. 
   RealType BoxSizes[MAX_N_OMEGARANGES][4];          // for each omegarange a box (window: left  right  bottom top) that defines the spots to include during indexing [micron]  
   RealType OmegaRanges[MAX_N_OMEGARANGES][2];       // Omegaranges: min, max [degrees], multiple possible.
   char OutputFolder[1024];        // output folder    
   int NoOfOmegaRanges;            // Automaticly set from Omegaranges (not explicit input by user)
   char SpotsFileName[1024];       // filename containing observed spots (see top for definition of columns) 
   char IDsFileName [1024];        // filename containing the spot-ids that will be used for indexing
   int UseFriedelPairs;            // 0=do not use friedelpairs  1=try to use friedelpairs
}; 


////////////////////////////////////////////////////////////////////////////////
int
ReadParams(
   char FileName[],
   struct TParams * Params)
   
{
   #define MAX_LINE_LENGTH 1024
            
   FILE *fp;
   char line[MAX_LINE_LENGTH];   
   char dummy[MAX_LINE_LENGTH];
   char *str;
   int NrOfBoxSizes = 0;
   int cmpres;   
   int NoRingNumbers = 0; // should be equal to Params->NrOfRings
   Params->NrOfRings = 0;
   Params->NoOfOmegaRanges = 0;
   
   fp = fopen(FileName, "r");
   if (fp==NULL) {
      printf("Cannot open file: %s.\n", FileName);
      return(1);
   }   
   
   // now get the params: format: "string" value(s)
   while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {

      str = "RingNumbers ";
      cmpres = strncmp(line, str, strlen(str));
      if (cmpres == 0) {
         sscanf(line, "%s %d", dummy, &(Params->RingNumbers[NoRingNumbers]) );
         NoRingNumbers++;
         continue;
      }   

      str = "CellStruct ";
      cmpres = strncmp(line, str, strlen(str));
      if (cmpres == 0) {
         sscanf(line, "%s %d", dummy, &(Params->CellStruct) );
         continue;
      }   
      
      str = "LatticeConstant ";
      cmpres = strncmp(line, str, strlen(str));
      if (cmpres == 0) {
         sscanf(line, "%s %lf", dummy, &(Params->LatticeConstant) );
         continue;
      }       
                             
      str = "Wavelength ";                           
      cmpres = strncmp(line, str, strlen(str));
      if (cmpres == 0) {
         sscanf(line, "%s %lf", dummy, &(Params->Wavelength) );
         continue;         
      }
      
      str = "Distance ";
      cmpres = strncmp(line, str, strlen(str));   
      if (cmpres == 0) {
         sscanf(line, "%s %lf", dummy, &(Params->Distance) );
         continue;         
      }     

      str = "Rsample ";
      cmpres = strncmp(line, str, strlen(str));      
      if ( cmpres == 0) {
         sscanf(line, "%s %lf", dummy, &(Params->Rsample) );
         continue;         
      }     
      
      str = "Hbeam ";
      cmpres = strncmp(line, str, strlen(str));             
      if ( cmpres == 0) {
         sscanf(line, "%s %lf", dummy, &(Params->Hbeam) );
         continue;         
      }

      str = "StepsizePos ";
      cmpres = strncmp(line, str, strlen(str));         
      if (cmpres == 0) {
         sscanf(line, "%s %lf", dummy, &(Params->StepsizePos) );
         continue;         
      }     
      
      str = "StepsizeOrient ";
      cmpres = strncmp(line, str, strlen(str));         
      if (cmpres == 0) {      
         sscanf(line, "%s %lf", dummy, &(Params->StepsizeOrient) );
         continue;         
      }     

      str = "MarginOme ";
      cmpres = strncmp(line, str, strlen(str));         
      if (cmpres == 0) {      
         sscanf(line, "%s %lf", dummy, &(Params->MarginOme) );
         continue;         
      }     

      str = "MarginRadius ";
      cmpres = strncmp(line, str , strlen(str));   
      if (cmpres == 0) {         
         sscanf(line, "%s %lf", dummy, &(Params->MarginRad) );
         continue;         
      }
      
      str = "MarginRadial ";
      cmpres = strncmp(line, str, strlen(str));         
      if (cmpres == 0) {         
         sscanf(line, "%s %lf", dummy, &(Params->MarginRadial) );
         continue;         
      }
      
      str = "EtaBinSize ";
      cmpres = strncmp(line, str, strlen(str));         
      if (cmpres == 0) {         
         sscanf(line, "%s %lf", dummy, &(Params->EtaBinSize) );
         continue;         
      }

      str = "OmeBinSize ";
      cmpres = strncmp(line, str, strlen(str));         
      if (cmpres == 0) {         
         sscanf(line, "%s %lf", dummy, &(Params->OmeBinSize) );
         continue;         
      }      

      str = "MinMatchesToAcceptFrac ";      
      cmpres = strncmp(line, str, strlen(str));            
      if (cmpres == 0) {   
         sscanf(line, "%s %lf", dummy, &(Params->MinMatchesToAcceptFrac) );
         continue;         
      }  
      
      str = "ExcludePoleAngle ";
      cmpres = strncmp(line, str, strlen(str));             
      if (cmpres == 0) {   
         sscanf(line, "%s %lf", dummy, &(Params->ExcludePoleAngle) );
         continue;         
      }  
            
      str = "RingRadii ";
      cmpres = strncmp(line, str, strlen(str));             
      if (cmpres == 0) {
         sscanf(line, "%s %lf", dummy, &(Params->RingRadiiUser[Params->NrOfRings]));  
         Params->NrOfRings = Params->NrOfRings + 1;           
         continue;
      }
      
      str = "OmegaRange ";
      cmpres = strncmp(line, str, strlen(str));             
      if (cmpres == 0) {       
         sscanf(line, "%s %lf %lf", dummy, &(Params->OmegaRanges[Params->NoOfOmegaRanges][0]),
                                         &(Params->OmegaRanges[Params->NoOfOmegaRanges][1]));
         (Params->NoOfOmegaRanges)++;
         continue;                                                            
      }
 
      str = "BoxSize ";
      cmpres = strncmp(line, str, strlen(str));             
      if (cmpres == 0) {       
         sscanf(line, "%s %lf %lf %lf %lf", dummy, &(Params->BoxSizes[NrOfBoxSizes][0]),
                                               &(Params->BoxSizes[NrOfBoxSizes][1]),
                                               &(Params->BoxSizes[NrOfBoxSizes][2]),
                                               &(Params->BoxSizes[NrOfBoxSizes][3]));
         NrOfBoxSizes++;
         continue;                                                        
      }
            
      str = "SpotsFileName ";
      cmpres = strncmp(line, str, strlen(str));             
      if (cmpres == 0) {            
         sscanf(line, "%s %s", dummy, Params->SpotsFileName );
         continue;         
      }
      
      str = "IDsFileName ";
      cmpres = strncmp(line, str, strlen(str));             
      if (cmpres == 0) {            
         sscanf(line, "%s %s", dummy, Params->IDsFileName  );
         continue;         
      }
      
      str = "MarginEta ";      
      cmpres = strncmp(line, str, strlen(str));            
      if (cmpres == 0) {   
         sscanf(line, "%s %lf", dummy, &(Params->MarginEta) );
         continue;         
      }

      str = "UseFriedelPairs ";      
      cmpres = strncmp(line, str, strlen(str));    
      if (cmpres == 0) {   
         sscanf(line, "%s %d", dummy, &(Params->UseFriedelPairs) );
         continue;         
      }                 
      
      str = "OutputFolder ";      
      cmpres = strncmp(line, str, strlen(str));    
      if (cmpres == 0) {   
         sscanf(line, "%s %s", dummy, Params->OutputFolder );
         continue;         
      }   
            
      // if string is empty 
      str = "";
      cmpres = strncmp(line, str, strlen(str));
      if (cmpres == 0) {
         continue;
      }
      
      // if string not recognized: print warning all other cases
      printf("Warning: skipping line in parameters file:\n");
      printf("%s\n", line);
   }
   
   // make a Params->RingRadii for internal use: ringno is directly the index in array (RingRadii[5] = ringradius from ring 5)
   int i;
   for (i = 0 ; i < MAX_N_RINGS ; i++ ) { Params->RingRadii[i] = 0; }
   for (i = 0 ; i < Params->NrOfRings ; i++ ) {
      Params->RingRadii[Params->RingNumbers[i]] = Params->RingRadiiUser[i]; 
   }
   
   return(0);
}


////////////////////////////////////////////////////////////////////////////////
int
WriteParams(
   char FileName[],
   struct TParams Params)
{
   
   FILE *fp;
   int i;
   
   fp = fopen(FileName, "w");
   if (fp==NULL) {
      printf("Cannot open file: %s.\n", FileName);
      return (1);
   }   
   
   int col1w = -22;
   
   // use alphabetic order
   for (i=0; i<Params.NoOfOmegaRanges; i++) {
      fprintf(fp, "%*s %lf %lf %lf %lf\n", col1w, "BoxSize ", Params.BoxSizes[i][0], Params.BoxSizes[i][1], Params.BoxSizes[i][2], Params.BoxSizes[i][3]);   
   }      
   
   fprintf(fp, "%*s %d\n", col1w, "CellStruct ", Params.CellStruct);   
   fprintf(fp, "%*s %lf\n", col1w, "Distance ", Params.Distance);
   fprintf(fp, "%*s %lf\n", col1w, "EtaBinSize ", Params.EtaBinSize);
   fprintf(fp, "%*s %lf\n", col1w, "ExcludePoleAngle ", Params.ExcludePoleAngle);
   fprintf(fp, "%*s %lf\n", col1w, "Hbeam ", Params.Hbeam);
   fprintf(fp, "%*s %s\n", col1w, "IDsFileName ", Params.IDsFileName );
   fprintf(fp, "%*s %lf\n", col1w, "LatticeConstant ", Params.LatticeConstant);
   fprintf(fp, "%*s %lf\n", col1w, "MarginEta ", Params.MarginEta);     
   fprintf(fp, "%*s %lf\n", col1w, "MarginOme ", Params.MarginOme);                   
   fprintf(fp, "%*s %lf\n", col1w, "MarginRadius ", Params.MarginRad);
   fprintf(fp, "%*s %lf\n", col1w, "MarginRadial ", Params.MarginRadial);
   fprintf(fp, "%*s %lf\n", col1w, "MinMatchesToAcceptFrac ", Params.MinMatchesToAcceptFrac);
   fprintf(fp, "%*s %lf\n", col1w, "OmeBinSize ", Params.OmeBinSize);
  

   for (i=0; i<Params.NoOfOmegaRanges; i++) {
      fprintf(fp, "%*s %lf %lf\n", col1w, "OmegaRange ", Params.OmegaRanges[i][0], Params.OmegaRanges[i][1]);   
   }

   fprintf(fp, "%*s %s\n", col1w, "OutputFolder ", Params.OutputFolder);
  
   for (i=0; i<Params.NrOfRings; i++) {
      fprintf(fp, "%*s %lf\n", col1w, "RingRadii ", Params.RingRadiiUser[i]);
   }
   
   for (i=0; i<Params.NrOfRings; i++) {
      fprintf(fp, "%*s %d\n", col1w, "RingNumbers ", Params.RingNumbers[i]);
   }   

   fprintf(fp, "%*s %lf\n", col1w, "Rsample ", Params.Rsample);
   fprintf(fp, "%*s %s\n", col1w, "SpotsFileName ", Params.SpotsFileName);
   fprintf(fp, "%*s %lf\n", col1w, "StepsizePos ", Params.StepsizePos);
   fprintf(fp, "%*s %lf\n", col1w, "StepsizeOrient ", Params.StepsizeOrient);
   fprintf(fp, "%*s %d\n", col1w, "UseFriedelPairs ",  Params.UseFriedelPairs);      
   fprintf(fp, "%*s %lf\n", col1w, "Wavelength ",  Params.Wavelength);   
      
   fclose(fp);
   
   return(0);
}


////////////////////////////////////////////////////////////////////////////////
void
CreateNumberedFilename(
   char stem[1000],
   int aNumber,
   char ext[10],
   char fn[1000+10+4]) 
{
   sprintf(fn, "%s%04d%s", stem, aNumber, ext);
}


////////////////////////////////////////////////////////////////////////////////
// create a filename with a number: ie file0003.txt
// the number of digits is given by numberOfDigits (in above example 4)
//
void
CreateNumberedFilenameW(
   char stem[1000],
   int aNumber,
   int numberOfDigits,
   char ext[10],
   char fn[1000+10+numberOfDigits+1]) 
{
   sprintf(fn, "%s%0*d%s", stem, numberOfDigits, aNumber, ext);
}




////////////////////////////////////////////////////////////////////////////////  
void   
GenerateRingInfoBCC(
   RealType LatticeConstant,
   RealType Wavelength,
   RealType RingTtheta[],
   int   RingMult[],
   int   RingHKL[][3]) 

{
   RealType theta;
   
   RingTtheta[0] = 0;
   
   CalcTheta(1, 1, 0, LatticeConstant, Wavelength, &theta);
   RingTtheta[1] = 2*theta;
   
   CalcTheta(2, 0, 0, LatticeConstant, Wavelength, &theta);
   RingTtheta[2] = 2*theta;
       
   CalcTheta(2, 1, 1, LatticeConstant, Wavelength, &theta);
   RingTtheta[3] = 2*theta;
       
   CalcTheta(2, 2, 0, LatticeConstant, Wavelength, &theta);
   RingTtheta[4] = 2*theta;
   
   CalcTheta(3, 1, 0, LatticeConstant, Wavelength, &theta);
   RingTtheta[5] = 2*theta; 
   
   RingMult[0] = 0;
   RingMult[1] = 12;
   RingMult[2] = 6;
   RingMult[3] = 24;
   RingMult[4] = 12;
   RingMult[5] = 24;         
   
   RingHKL[0][0] = 0;  
   RingHKL[0][1] = 0; 
   RingHKL[0][2] = 0; 
   
   RingHKL[1][0] = 1;  
   RingHKL[1][1] = 1; 
   RingHKL[1][2] = 0;

   RingHKL[2][0] = 2;  
   RingHKL[2][1] = 0; 
   RingHKL[2][2] = 0;
   
   RingHKL[3][0] = 2;  
   RingHKL[3][1] = 1; 
   RingHKL[3][2] = 1;
   
   RingHKL[4][0] = 2;  
   RingHKL[4][1] = 2; 
   RingHKL[4][2] = 0;
   
   RingHKL[5][0] = 3;  
   RingHKL[5][1] = 1; 
   RingHKL[5][2] = 0;                 

}

 
////////////////////////////////////////////////////////////////////////////////
void
PrintRingInfo(
   RealType RingTtheta[],
   int RingMult[],
   int RingHKL[][3],
   int RingNos[],
   int nRings) 
   
{
   int i, RingNo;
   
   printf("Ring info:\n");
   printf("RingNo  h k l  mult  2-theta \n");  
   for (i=0 ; i < nRings ; i++){
      RingNo = RingNos[i];
      printf("%6d  %d %d %d   %2d   %lf\n", RingNo, RingHKL[RingNo][0], RingHKL[RingNo][1], RingHKL[RingNo][2], RingMult[RingNo], RingTtheta[RingNo] );
   }
   printf("\n");
}    


////////////////////////////////////////////////////////////////////////////////  
void   
GenerateRingInfoFCC(
   RealType LatticeConstant,
   RealType Wavelength,
   RealType RingTtheta[],
   int   RingMult[],
   int   RingHKL[][3]) 

{
   RealType theta;
   
   RingTtheta[0] = 0;
   CalcTheta(1, 1, 1, LatticeConstant, Wavelength, &theta);
   RingTtheta[1] = 2*theta;
   CalcTheta(2, 0, 0, LatticeConstant, Wavelength, &theta);
   RingTtheta[2] = 2*theta;    
   CalcTheta(2, 2, 0, LatticeConstant, Wavelength, &theta);
   RingTtheta[3] = 2*theta;    
   CalcTheta(3, 1, 1, LatticeConstant, Wavelength, &theta);
   RingTtheta[4] = 2*theta;
   CalcTheta(2, 2, 2, LatticeConstant, Wavelength, &theta);
   RingTtheta[5] = 2*theta; 
   
   RingMult[0] = 0;
   RingMult[1] = 8;
   RingMult[2] = 6;
   RingMult[3] = 12;
   RingMult[4] = 24;
   RingMult[5] = 8;         
   
   RingHKL[0][0] = 0;  
   RingHKL[0][1] = 0; 
   RingHKL[0][2] = 0; 
   
   RingHKL[1][0] = 1;  
   RingHKL[1][1] = 1; 
   RingHKL[1][2] = 1;

   RingHKL[2][0] = 2;  
   RingHKL[2][1] = 0; 
   RingHKL[2][2] = 0;
   
   RingHKL[3][0] = 2;  
   RingHKL[3][1] = 2; 
   RingHKL[3][2] = 0;
   
   RingHKL[4][0] = 3;  
   RingHKL[4][1] = 1; 
   RingHKL[4][2] = 1;             

   RingHKL[5][0] = 2;  
   RingHKL[5][1] = 2; 
   RingHKL[5][2] = 2;

}


////////////////////////////////////////////////////////////////////////////////
void   
GenerateRingInfo(
   int cellstruct,  // 1=bcc 2=fcc
   RealType LatticeConstant,
   RealType Wavelength,
   RealType RingTtheta[],
   int   RingMult[],
   int   RingHKL[][3]) 
{
    switch  ( cellstruct )
    {
       case 1 : GenerateRingInfoBCC(LatticeConstant, Wavelength, RingTtheta, RingMult, RingHKL);
                break;
       case 2 : GenerateRingInfoFCC(LatticeConstant, Wavelength, RingTtheta, RingMult, RingHKL);
                break;      
    }
}


////////////////////////////////////////////////////////////////////////////////
RealType
CalcAvgIA(RealType *Arr, int n)
{
   RealType total = 0;
   int nnum = 0;
   int i;
   
   for (i=0 ; i < n ; i++) {
       if (Arr[i] == 999)     continue;  // skip special number
       total = total + fabs(Arr[i]);
       nnum++;
   }
   
   if (nnum == 0)  
      return 0;
   else 
      return total / nnum;

}

////////////////////////////////////////////////////////////////////////////////
// selects the best grain (smallest avg IA) and writes it to a file, including
// all spots belonging to this grain.
//
int
WriteBestMatch(
   char *FileName,
   RealType **GrainMatches, 
   int ngrains,
   RealType **AllGrainSpots,  // N_COL_GRAINSPOTS 
   int nrows,
   char *FileName2)
   
{
   int r, g, c;
   RealType smallestIA = 99999;
   int bestGrainIdx = -1;
   
   // find the smallest ia
   for ( g = 0 ; g < ngrains ; g++ ) {
      if ( GrainMatches[g][15] < smallestIA ) {
         
         smallestIA = GrainMatches[g][15];
         bestGrainIdx = g;
      }
   }
   

   // if found, write to file
   if (bestGrainIdx != -1) {
      FILE *fp;
      FILE *fp2;
  
      fp = fopen(FileName, "w");
      fp2 = fopen(FileName2,"w");
      if (fp2==NULL) {
           printf("Cannot open file: %s\n", FileName2);
           return(1);
      }
      if (fp==NULL) {
         printf("Cannot open file: %s\n", FileName);
         return(1);
      }   
  
      // write grain info (orient and pos)
      fprintf(fp, "o11 o12 o13 o21 o22 o23 o31 o32 o33 x y z NTheor NMatched MatchId avgIA\n");      
      for (c = 0; c < N_COL_GRAINMATCHES; c++) {
         fprintf(fp, "%lf ", GrainMatches[bestGrainIdx][c]);
      }
      fprintf(fp, "\n\n");
      // write spots
      RealType bestGrainID =  GrainMatches[bestGrainIdx][14];
      fprintf(fp2, "%lf\n%lf\n",bestGrainID,bestGrainID);
       fprintf(fp2,"%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", GrainMatches[bestGrainIdx][15], GrainMatches[bestGrainIdx][0], GrainMatches[bestGrainIdx][1], GrainMatches[bestGrainIdx][2], GrainMatches[bestGrainIdx][3], GrainMatches[bestGrainIdx][4], GrainMatches[bestGrainIdx][5], GrainMatches[bestGrainIdx][6], GrainMatches[bestGrainIdx][7], GrainMatches[bestGrainIdx][8], GrainMatches[bestGrainIdx][9], GrainMatches[bestGrainIdx][10], GrainMatches[bestGrainIdx][11], GrainMatches[bestGrainIdx][12], GrainMatches[bestGrainIdx][13]);
      fprintf(fp, "No Dummy YCalc YObs YDiff ZCalc ZObs ZDiff WCalc WObs WDiff RadRef RadObs RadDiff SpotId MatchNo IA\n");
      for (r = 0 ; r < nrows ; r++ ) {
         if (AllGrainSpots[r][15] == bestGrainID ) {
            for (c = 0; c < N_COL_GRAINSPOTS; c++) {
               fprintf(fp, "%14lf ", AllGrainSpots[r][c]);
                if (c!=1) {
                    fprintf(fp2,"%14lf, ", AllGrainSpots[r][c]);
                }
            }
            fprintf(fp, "\n");
            fprintf(fp2, "\n");
         }
      }
      fclose(fp);
      fclose(fp2);
   }
   else {
       FILE *fp2;
       fp2 = fopen(FileName2,"w");
       fclose(fp2);
   }
   return (0);
}



////////////////////////////////////////////////////////////////////////////////
void
CalcIA(
   RealType **GrainMatches, //N_COL_GRAINMATCHES,
   int ngrains,
   RealType **AllGrainSpots, // N_COL_GRAINSPOTS
   RealType distance)
   
   
{
   RealType *IAgrainspots;
   int r, g;
   RealType g1x, g1y, g1z;
   RealType x1, y1, z1, w1, x2, y2, z2, w2, gv1x, gv1y, gv1z, gv2x, gv2y, gv2z;
   
   // Calc for AllGrainSpots Internal angle for each spot
   // nr dummy yo yt yd z... w... rad (12) spotid matchnr
   // calc ia between the 2 g-vectors 
   // yzw1 is the vector to observed spot
   // yzw2 is the vector to theoretical spot
   // substract grain origin to get a vector from the grain origin (instead of 
   // grain centre)
   
   int nspots;    
   int rt = 0;    // rt = current row in all grain spots matrix
   
  // for faster calc, alloc only once (NB should be large enough to hold all spots of 1 grain) 
  IAgrainspots = malloc(10000 * sizeof(* IAgrainspots));   
   for ( g = 0 ; g < ngrains ; g++ ) {
      nspots = GrainMatches[g][12]; // no of spots stored (also the theoretical spots not found, just skip those)
      for (r=0 ; r < nspots ; r++) {
        
         // skip negative ids
         if (AllGrainSpots[rt][0] < 0) {
            AllGrainSpots[rt][16] = 999;   // a special number (avoid 0)
            IAgrainspots[r] = AllGrainSpots[rt][16];  
            rt++; 
            continue;
         } 
         
         x1 = distance;
         x2 = distance; 
         
         y1 = AllGrainSpots[rt][2];
         y2 = AllGrainSpots[rt][3];
         
         z1 = AllGrainSpots[rt][5];
         z2 = AllGrainSpots[rt][6];
         
         w1 = AllGrainSpots[rt][8];
         w2 = AllGrainSpots[rt][9];
         
         // grain pos
         g1x = GrainMatches[g][9];
         g1y = GrainMatches[g][10];
         g1z = GrainMatches[g][11];
         
         spot_to_gv_pos(x1, y1, z1, w1, g1x, g1y, g1z, &gv1x, &gv1y, &gv1z);
         spot_to_gv_pos(x2, y2, z2, w2, g1x, g1y, g1z, &gv2x, &gv2y, &gv2z);      
         
         CalcInternalAngle(gv1x, gv1y, gv1z, gv2x, gv2y, gv2z, &AllGrainSpots[rt][16]);
         IAgrainspots[r] = AllGrainSpots[rt][16];
         rt++; 
      }
      GrainMatches[g][15] = CalcAvgIA(IAgrainspots, nspots);      
   }
   free(IAgrainspots);   
}


////////////////////////////////////////////////////////////////////////////////
void MakeFullFileName(char* fullFileName, char* aPath, char* aFileName)
{
   if (aPath[0] == '\0' )  {
       strcpy(fullFileName, aFileName); 
   }
   else {
      // else concat path and fn
      strcpy(fullFileName, aPath);
      strcat(fullFileName, "/");
      strcat(fullFileName, aFileName);
   }
}






////////////////////////////////////////////////////////////////////////////////
// tries to find grains (orientations for each spot in spotIDs) 
int 
DoIndexing(
//     ObsSpotsLab, global
//   int nspots,
   int SpotIDs,
   int nSpotIDs,
   struct TParams Params )

{
   time_t start, end;
   double dif;  
      
   RealType HalfBeam = Params.Hbeam /2 ;      
   RealType MinMatchesToAccept;  
   RealType ga, gb, gc;  
   int   nTspots;
   int   bestnMatchesIsp, bestnMatchesRot, bestnMatchesPos;
   int   bestnTspotsIsp, bestnTspotsRot, bestnTspotsPos;
   int   matchNr;
   int   nOrient;
   RealType hklnormal[3]; 
   RealType Displ_y;
   RealType Displ_z;
   int   or;
   int   sp;
   int   nMatches;
   int   r,c, i;
   RealType y0_vector[MAX_N_STEPS];
   RealType z0_vector[MAX_N_STEPS];
   int   nPlaneNormals;
   int   hkl[3];
   RealType g1, g2, g3;
   int   isp;
   RealType xi, yi, zi;
   int   n_max, n_min, n;
   RealType y0, z0; 
   RealType RingTtheta[MAX_N_RINGS];
   int   RingMult[MAX_N_RINGS];
   int   RingHKL[MAX_N_RINGS][3];
   int   orDelta, ispDelta, nDelta;
   RealType fracMatches;
   int   rownr;
   int   SpotRowNo;
   int usingFriedelPair; // 0 = not, 1 = yes, using fridelpair
   RealType **BestMatches;
   
   RealType omemargins[181];
   RealType etamargins[MAX_N_RINGS];
   char fn[1000]; 
   char ffn[1000]; // full fn
    char fn2[1000];
   char ffn2[1000];
   
   // output matrix: grains found for 1 spot
//   RealType GrainMatches[MAX_N_MATCHES][N_COL_GRAINMATCHES];    // orientat (9), pos (3), nTheorectical, nMatched, machtnr
   RealType **GrainMatches; //[MAX_N_MATCHES][N_COL_GRAINMATCHES];    // orientat (9), pos (3), nTheorectical, nMatched, machtnr   
   RealType **TheorSpots;   // theoretical spots (2d matrix, see top for columns)   
   RealType **GrainSpots;   // spots found for 1 grain ( pos + orient )     
   RealType **AllGrainSpots; // for each match, the corresponding spots and their differences:
                         // nr dummy yo yt yd z... w... rad (12) spotid matchnr

   // allocate output matrices
   
   // grainspots for all grains (from 1 starting spot)
   int nRowsOutput = MAX_N_MATCHES * 2 * n_hkls;
   AllGrainSpots = allocMatrix(nRowsOutput, N_COL_GRAINSPOTS);
   if (AllGrainSpots == NULL ) {
      printf("Memory error: could not allocate memory for output matrix. Memory full?\n");
      return 1;
   }      
      
   // grain matches (pos + orient)
   GrainMatches = allocMatrix(MAX_N_MATCHES, N_COL_GRAINMATCHES); 
   if (GrainMatches == NULL ) {
      printf("Memory error: could not allocate memory for output matrix. Memory full?\n");
      return 1;
   }      
   
   // grainspots for 1 grain
   int nRowsPerGrain = 2 * n_hkls;
   GrainSpots = allocMatrix(nRowsPerGrain, N_COL_GRAINSPOTS ); 
   

   // theorspots
   TheorSpots = allocMatrix(nRowsPerGrain, N_COL_THEORSPOTS);
   if (TheorSpots == NULL ) {
      printf("Memory error: could not allocate memory for output matrix. Memory full?\n");
      return 1;
   }      
   
   // bestmatches
   BestMatches = allocMatrix(nSpotIDs, 5);
   if (BestMatches == NULL ) {
      printf("Memory error: could not allocate memory for output matrix. Memory full?\n");
      return 1;
   }      
   
   // store information for each ring
   GenerateRingInfo(Params.CellStruct, Params.LatticeConstant, Params.Wavelength, RingTtheta, RingMult, RingHKL);
   PrintRingInfo(RingTtheta, RingMult, RingHKL, Params.RingNumbers, Params.NrOfRings );
   
   // pre calc omega margins
   for ( i = 1 ; i < 180 ; i++) {
      omemargins[i] = Params.MarginOme + ( 0.5 * Params.StepsizeOrient / fabs(sin(i * deg2rad))); 
   }
   omemargins[0] = omemargins[1]; // officially undefined
   omemargins[180] = omemargins[1];   
    
   // etamargins 
   for ( i = 0 ; i < MAX_N_RINGS ; i++) {
      if ( Params.RingRadii[i] == 0)  { 
         etamargins[i] = 0; 
      }
      else {
         etamargins[i] = rad2deg * atan(Params.MarginEta/Params.RingRadii[i]) + 0.5 * Params.StepsizeOrient;
      }   
   }
   
//   WriteArrayF("omemargins.txt", omemargins, 181);
//   WriteArrayF("etamargins.txt", etamargins, Params.NrOfRings+1);   
        
   int SpotIDIdx;
   printf("Starting indexing...\n");
   for ( SpotIDIdx=0 ; SpotIDIdx < nSpotIDs ; SpotIDIdx++) {
//   for ( SpotIDIdx=1 ; SpotIDIdx < 2 ; SpotIDIdx++) {   
      time (&start);   
      matchNr = 0;   // grain match no
      rownr = 0; // row no in the output file of all the spots of the grains found of this spot
  
      // find row for spot to check
      RealType SpotID = SpotIDs;
      FindInMatrix(&ObsSpotsLab[0][0], n_spots, N_COL_OBSSPOTS, 4, SpotID, &SpotRowNo);
    
      if (SpotRowNo == -1) {
         printf("WARNING: SpotId %lf not found in spots file! Ignoring this spotID.\n", SpotID);
         continue;
      }
      
      RealType ys     = ObsSpotsLab[SpotRowNo][0];   
      RealType zs     = ObsSpotsLab[SpotRowNo][1];
      RealType omega  = ObsSpotsLab[SpotRowNo][2];  
      RealType RefRad = ObsSpotsLab[SpotRowNo][3];
      RealType eta    = ObsSpotsLab[SpotRowNo][6];
      RealType ttheta = ObsSpotsLab[SpotRowNo][7];
      int   ringnr = (int) ObsSpotsLab[SpotRowNo][5];  
      
      hkl[0] = RingHKL[ringnr][0];
      hkl[1] = RingHKL[ringnr][1];      
      hkl[2] = RingHKL[ringnr][2];      
      

      printf("\n--------------------------------------------------------------------------\n");
      printf("Spot number being processed %d of %d\n\n",SpotIDIdx, nSpotIDs);
      printf("%8s %10s %9s %9s %9s %9s %9s %7s\n", "SpotID", "SpotRowNo", "ys", "zs", "omega", "eta", "ttheta", "ringno");
      printf("%8.0f %10d %9.2f %9.2f %9.3f %9.3f %9.3f %7d\n\n", SpotID, SpotRowNo, ys, zs, omega, eta, ttheta, ringnr);
      
      // Generate 'ideal spots on the ring', NB ttheta input is ttheta of ideal ring
      nPlaneNormals = 0;
      if (Params.UseFriedelPairs == 1) {
         usingFriedelPair = 1;  // 1= true 0 = false
         GenerateIdealSpotsFriedel(ys, zs, RingTtheta[ringnr], eta, omega, ringnr, 
              Params.RingRadii[ringnr], Params.Rsample, Params.Hbeam, Params.MarginOme, Params.MarginRadial, 
              y0_vector, z0_vector, &nPlaneNormals);
         
         // check of friedelpair was found (nplanenormals > 0). If not try to find the mixed friedelpair     
         if (nPlaneNormals == 0 ) {
            GenerateIdealSpotsFriedelMixed(ys, zs, RingTtheta[ringnr], eta, omega, ringnr, 
              Params.RingRadii[ringnr], Params.Distance, Params.Rsample, Params.Hbeam, Params.StepsizePos, 
              Params.MarginOme, Params.MarginRadial, Params.MarginEta,
              y0_vector, z0_vector, &nPlaneNormals);
         }
      }

      // if friedelpair was not found (nPlaneNormals = 0): do full search of plane      
      if ( nPlaneNormals == 0 ) {
         usingFriedelPair = 0; // affects skipping
         printf("Trying all plane normals.\n");
         GenerateIdealSpots(ys, zs, RingTtheta[ringnr], eta, Params.RingRadii[ringnr], Params.Rsample, Params.Hbeam, Params.StepsizePos, y0_vector, z0_vector, &nPlaneNormals);
      }
      
      printf("No of Plane normals: %d\n\n", nPlaneNormals);
//      WriteArrayF("PlaneNY0.txt", y0_vector, nPlaneNormals);
//      WriteArrayF("PlaneNZ0.txt", z0_vector, nPlaneNormals);
      
      // try each plane normal
      bestnMatchesIsp = -1;
      bestnTspotsIsp = 0;
      isp = 0;
      while (isp < nPlaneNormals) {
         y0 = y0_vector[isp];
         z0 = z0_vector[isp];
         MakeUnitLength(Params.Distance, y0, z0, &xi, &yi, &zi );         
         spot_to_gv(xi, yi, zi, omega,  &g1, &g2, &g3);
         
         hklnormal[0] = g1;
         hklnormal[1] = g2;  
         hklnormal[2] = g3;
         
         // Generate candidate orientations for this spot (= plane)
         GenerateCandidateOrientationsF(hkl, hklnormal, Params.StepsizeOrient, OrMat, &nOrient);
        
         // try each rotation
         bestnMatchesRot = -1;
         bestnTspotsRot = 0;         
         or = 0;
         orDelta = 1; 
         while (or < nOrient) {
            // convert from row to matrix format
//            for (i = 0 ;  i < 9 ; i ++) { OrientMatrix[i/3][i%3] =  OrMat[or][i]; }
         
            CalcDiffrSpots_Furnace(OrMat[or], Params.LatticeConstant, Params.Wavelength , Params.Distance, Params.RingRadii, Params.OmegaRanges, Params.BoxSizes, Params.NoOfOmegaRanges, Params.ExcludePoleAngle, TheorSpots, &nTspots);
#ifdef DEBUG            
            printf("nTspots: %d\n", nTspots);
#endif            
            MinMatchesToAccept = nTspots * Params.MinMatchesToAcceptFrac;
            
            // step in sample
            bestnMatchesPos = -1;      
            bestnTspotsPos =  0;      
            calc_n_max_min(xi, yi, ys, y0, Params.Rsample, Params.StepsizePos, &n_max, &n_min);
#ifdef DEBUG            
            printf("n_min n_max: %d %d\n", n_min, n_max);   
#endif
            n = n_min;
            
            // step in the sample
            while (n <= n_max) {
               spot_to_unrotated_coordinates(xi, yi, zi, ys, zs, y0, z0, Params.StepsizePos, n, omega, &ga, &gb, &gc );
               if (fabs(gc) > HalfBeam) {
                   n++; 
                   continue; // outside sample               
               }   
 
               // displace spots due to diplacement of grain in sample
               for (sp = 0 ; sp < nTspots ; sp++) {
                  displacement_spot_needed_COM(ga, gb, gc, TheorSpots[sp][3], TheorSpots[sp][4], 
                      TheorSpots[sp][5], TheorSpots[sp][6], &Displ_y, &Displ_z );    
                                                          
                  TheorSpots[sp][10] = TheorSpots[sp][4] +  Displ_y;
                  TheorSpots[sp][11] = TheorSpots[sp][5] +  Displ_z;
                  CalcEtaAngle( TheorSpots[sp][10], TheorSpots[sp][11], &TheorSpots[sp][12] ); // correct eta for displaced spot
                  TheorSpots[sp][13] = sqrt(TheorSpots[sp][10] * TheorSpots[sp][10] + TheorSpots[sp][11] * TheorSpots[sp][11]) - 
                               Params.RingRadii[(int)TheorSpots[sp][9]];  // new eta, due to displ spot                                   
               }
         
               // compare theoretical spots with experimental spots
               CompareSpots(TheorSpots, nTspots, ObsSpotsLab, RefRad, 
                        Params.MarginRad, Params.MarginRadial, etamargins, omemargins, 
                        &nMatches, GrainSpots);

               if (nMatches > bestnMatchesPos) {  
                  bestnMatchesPos = nMatches;
                  bestnTspotsPos = nTspots;   
               }

               // save output (// save the best match)
               if ( (nMatches > 0) &&
                    (matchNr < MAX_N_MATCHES) &&
                    (nMatches >= MinMatchesToAccept) ) { 
                  for (i = 0 ;  i < 9 ; i ++) { GrainMatches[matchNr][i] = OrMat[or][i/3][i%3]; }
               
                  GrainMatches[matchNr][9]  = ga;
                  GrainMatches[matchNr][10] = gb;
                  GrainMatches[matchNr][11] = gc;        
                  GrainMatches[matchNr][12] = nTspots;
                  GrainMatches[matchNr][13] = nMatches;
                  GrainMatches[matchNr][14] = matchNr+1;   
                   
                  // save difference for each spot and ID
                  for (r = 0 ; r < nTspots ; r++) {
                     for (c = 0 ; c < 15 ; c++) {
                        AllGrainSpots[rownr][c] = GrainSpots[r][c];
                     }
                     AllGrainSpots[rownr][15] = matchNr+1; // avoid 0 as matchnr
                     rownr++;                         
                  }
                  matchNr++;
               }
               
               // optimization: if in previous run the orientation gave low no of matches skip a few rotations
               nDelta = 1;  // default
               if (nTspots != 0) { 
                  fracMatches = (RealType)nMatches/nTspots;
                  if (fracMatches < 0.5) { nDelta = 10 - round(fracMatches * (10-1) / 0.5); } 
               }
                  
               n = n + nDelta;               
            }
            
            // save the best
            if (bestnMatchesPos > bestnMatchesRot) {  
               bestnMatchesRot = bestnMatchesPos;
               bestnTspotsRot = bestnTspotsPos;   
            }
            
#ifdef DEBUG     
            printf("isp nisp or nor npos nTheor nMatches: %d %d %d %d %d %d %d\n", isp, nPlaneNormals, or, nOrient, 2*n_max+1, bestnTspotsPos,  bestnMatchesPos );
#endif                     
            or = or + orDelta;            
         }
         
         // save the best               
         if (bestnMatchesRot > bestnMatchesIsp) {  
            bestnMatchesIsp = bestnMatchesRot;
            bestnTspotsIsp = bestnTspotsRot;   
         }             
         
         // optimization: if in previous run the plane normal gave low no of matches skip a few planenormals
         // if using Friedelpair: dont skip!
         ispDelta = 1;  // default
         if ((!usingFriedelPair) && (bestnTspotsRot != 0)) { 
            fracMatches = (RealType) bestnMatchesRot/bestnTspotsRot;
            if (fracMatches < 0.5) { ispDelta = 5 - round(fracMatches * (5-1) / 0.5); } // ispdelta between 1 and 5
         }

         printf("==> planenormal #pns #or #pos #Theor #Matches: %d %d %d %d %d %d\n", isp, nPlaneNormals, nOrient, 2*n_max+1, bestnTspotsRot, bestnMatchesRot);
         isp = isp + ispDelta;
         
         if (matchNr >= MAX_N_MATCHES) {
             printf("Warning: the number of grain matches exceeds maximum (%d). Not all output is saved!\n", MAX_N_MATCHES );
         }         
      }  // isp loop
      
      fracMatches = (RealType) bestnMatchesIsp/bestnTspotsIsp;      
      printf("\n==> Best Match: No_of_theoretical_spots No_of_spots_found fraction: %d %d %0.2f\n", bestnTspotsIsp, bestnMatchesIsp, fracMatches );
      BestMatches[SpotIDIdx][0] = SpotIDIdx+1;
      BestMatches[SpotIDIdx][1] = SpotID;      
      BestMatches[SpotIDIdx][2] = bestnTspotsIsp;      
      BestMatches[SpotIDIdx][3] = bestnMatchesIsp;      
      BestMatches[SpotIDIdx][4] = fracMatches;      
      
      time (&end);
      dif = difftime (end,start);
      printf("Time elapsed [s] [min]: %.0f %.0f\n", dif, dif/60);             

      // Calc IA for each spot and the avg per grain
      CalcIA(GrainMatches, matchNr, AllGrainSpots, Params.Distance );
       
      //CreateNumberedFilenameW("GrainSpots_", (int) SpotID, 9, ".txt", fn);
      //MakeFullFileName(ffn, Params.OutputFolder, fn);
      //WriteMatrixWithHeaderFp(ffn, AllGrainSpots, rownr, N_COL_GRAINSPOTS,
      //  "No Dummy YCalc YObs YDiff ZCalc ZObs ZDiff WCalc WObs WDiff RadRef RadObs RadDiff SpotId MatchNo IA");
      
      CreateNumberedFilenameW("BestGrain_", (int) SpotID, 9, ".txt", fn);
      MakeFullFileName(ffn, Params.OutputFolder, fn);
      CreateNumberedFilenameW("BestPos_", (int) SpotID, 9, ".csv", fn2);
      MakeFullFileName(ffn2, Params.OutputFolder, fn2);
      WriteBestMatch(ffn, GrainMatches, matchNr, AllGrainSpots, rownr, ffn2);
   }  // to next spot
   
   // write bestmatches
   //MakeFullFileName(ffn, Params.OutputFolder, "SpotScores.txt");
   //WriteMatrixWithHeaderFp(ffn, BestMatches, nSpotIDs, 5, "No SpotID NoSpotsTheor NoSpotsFound Frac");
   
   // free mem 
   FreeMemMatrix( GrainMatches, MAX_N_MATCHES);   
   FreeMemMatrix( TheorSpots, nRowsPerGrain);
   FreeMemMatrix( GrainSpots, nRowsPerGrain);   
   FreeMemMatrix( AllGrainSpots, nRowsOutput);
   FreeMemMatrix( BestMatches, nSpotIDs);
   
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
void
ConcatStr(char *str1, char *str2, char *resStr){
   strcpy(resStr, str1);
   strcat(resStr, str2);
}


////////////////////////////////////////////////////////////////////////////////
void
CalcDistanceIdealRing(RealType ObsSpotsLab[MAX_N_SPOTS][N_COL_OBSSPOTS], int nspots, RealType RingRadii[]  )
{
    int i;
    for (i = 0 ; i < nspots ; ++i)
    {
       RealType y = ObsSpotsLab[i][0];
       RealType z = ObsSpotsLab[i][1];
       RealType rad = sqrt(y*y + z*z);
       int ringno = (int) ObsSpotsLab[i][5];
       ObsSpotsLab[i][8] = rad - RingRadii[ringno]; 
    }

}


////////////////////////////////////////////////////////////////////////////////
//
// Show some statistics: how much spots in each bin.
void 
BinStats(int n_ring_bins, int n_eta_bins, int n_ome_bins)
{
   #define N_MAX_SPOT_IN_BIN 10

   int iRing, iEta, iOmega, i;
   int histogram[N_MAX_SPOT_IN_BIN+1];  // for 0
   
   for ( i = 0 ; i < N_MAX_SPOT_IN_BIN+1 ; i++) {
       histogram[i] = 0;
   }
      
   for ( iRing = 0 ; iRing < n_ring_bins ; iRing++ ){
       for ( iEta = 0 ; iEta < n_eta_bins ; iEta++) {
           for ( iOmega = 0 ; iOmega < n_ome_bins ; iOmega++) {
               int nspots = ndata[iRing][iEta][iOmega];
               if (nspots >  N_MAX_SPOT_IN_BIN) nspots = N_MAX_SPOT_IN_BIN;
               (histogram[nspots])++;    
           }
       }
   }
   
   printf("Bin statistics:\n");
   printf("%15s %13s\n", "#Spots_in_bin", "No_of_bins");
   for ( i = 0 ; i < N_MAX_SPOT_IN_BIN+1 ; i++) {
      printf("%15d %13d\n", i, histogram[i]);
   }
}


////////////////////////////////////////////////////////////////////////////////
int
MakeBins(
   RealType ObsSpots[][N_COL_OBSSPOTS],
   int nSpots,
   RealType omemargin0,
   RealType etamargin0,  // in micron!   
   RealType rotationstep,   
   RealType RingRadii[],
//   int NoOfRings,    // NB this is NOT the length of the array RingRadii. RingRadii has size: MAX_N_SPOTS (with a lot of zeros, for ringno that are not used)
   RealType etabinsize,
   RealType omebinsize)

{
   int *newarray;
   int *oldarray;
   int iEta, iOme, iEta0, iOme0;
   int rowno;
   
   // set global vars
   EtaBinSize = etabinsize;
   OmeBinSize = omebinsize;
   
   // get max ring no
   int HighestRingNo = 0;
   int i;
   for (i = 0 ; i < MAX_N_RINGS ; i++ ) { 
      if ( RingRadii[i] != 0) HighestRingNo = i;
   }
   // calculate no of bins (NB global vars)
   n_ring_bins = HighestRingNo;  
   n_eta_bins = ceil(360.0 / etabinsize);
   n_ome_bins = ceil(360.0 / omebinsize);   

   // init data
   int i1, i2, i3;
   data = malloc(n_ring_bins * sizeof(data));
   if (data == NULL ) {
      printf("Memory error: memory full?\n");
      return 1;
   }    
   for (i1 = 0 ; i1 < n_ring_bins ; i1++) {
       data[i1] = malloc(n_eta_bins * sizeof(data[i1]));
       if (data[i1] == NULL ) {
          printf("Memory error: memory full?\n");
          return 1;
       }
       for (i2 = 0 ; i2 < n_eta_bins ; i2++) {
           data[i1][i2] = malloc(n_ome_bins * sizeof(data[i1][i2]));
           if (data[i1][i2] == NULL ) {
              printf("Memory error: memory full?\n");
              return 1;
           }           
           for (i3 = 0 ; i3 < n_ome_bins ; i3++) {
               data[i1][i2][i3] = NULL;
           }
       } 
   }
  
   // init ndata
   ndata = malloc(n_ring_bins * sizeof(ndata));
   if (ndata == NULL ) {
      printf("Memory error: memory full?\n");
      return 1;
   }                                                              
   for (i1 = 0 ; i1 < n_ring_bins ; i1++) {
       ndata[i1] = malloc(n_eta_bins * sizeof(ndata[i1]));
       if (ndata[i1] == NULL ) {
          printf("Memory error: memory full?\n");
          return 1;
       }       
       for (i2 = 0 ; i2 < n_eta_bins ; i2++) {
           ndata[i1][i2] = malloc(n_ome_bins * sizeof(ndata[i1][i2]));
           if (ndata[i1][i2] == NULL ) {
              printf("Memory error: memory full?\n");
              return 1;
           }  
           for (i3 = 0 ; i3 < n_ome_bins ; i3++) {
               ndata[i1][i2][i3] = 0;
           }
       } 
   }

   // init maxdata
   maxndata = malloc(n_ring_bins * sizeof(maxndata));
   if (maxndata == NULL ) {
      printf("Memory error: memory full?\n");
      return 1;
   }    
   for (i1 = 0 ; i1 < n_ring_bins ; i1++) {
       maxndata[i1] = malloc(n_eta_bins * sizeof(maxndata[i1]));
       if (maxndata[i1] == NULL ) {
          printf("Memory error: memory full?\n");
          return 1;
       }       
       for (i2 = 0 ; i2 < n_eta_bins ; i2++) {
           maxndata[i1][i2] = malloc(n_ome_bins * sizeof(maxndata[i1][i2]));
           if (maxndata[i1][i2] == NULL ) {
              printf("Memory error: memory full?\n");
              return 1;
           }           
           for (i3 = 0 ; i3 < n_ome_bins ; i3++) {
               maxndata[i1][i2][i3] = 0;
           }
       } 
   }
   
   // put all spots in bins: a spot is also put in adjacent bins (depending on omega and eta ranges)
   for (rowno = 0 ; rowno < nSpots ; rowno++ ) {  
      int ringnr = (int) ObsSpots[rowno][5];
      RealType eta = ObsSpots[rowno][6];
      RealType omega = ObsSpots[rowno][2];

      int iRing = ringnr-1;
      if ( (iRing < 0) || (iRing > n_ring_bins-1) ) continue;
      if ( RingRadii[ringnr] == 0 ) continue;   // ring is not used

      RealType omemargin = omemargin0 + ( 0.5 * rotationstep / fabs(sin(eta * deg2rad)));
      RealType omemin = 180 + omega - omemargin;
      RealType omemax = 180 + omega + omemargin;
      int iOmeMin = floor(omemin / omebinsize);
      int iOmeMax = floor(omemax / omebinsize);

      RealType etamargin = rad2deg * atan(etamargin0/RingRadii[ringnr]) + 0.5 * rotationstep; 
      RealType etamin = 180 + eta - etamargin;
      RealType etamax = 180 + eta + etamargin;
      int iEtaMin = floor(etamin / etabinsize);     // iEta can be negative here , is corrected later       
      int iEtaMax = floor(etamax / etabinsize);

      // put this spot in all bins between mineta..maxeta and minomega..maxomega
      for ( iEta0 = iEtaMin ; iEta0 <= iEtaMax ; iEta0++) {
          // make sure iEta between 0 and n_bins_eta-1 (same for iOme)      
          iEta = iEta0 % n_eta_bins;
          if ( iEta < 0 ) iEta = iEta + n_eta_bins;  // rem can be neg!
          //printf("%d %d %d %f %d %f\n", iEta, iEtaMin, iEtaMax, etamargin, ringnr, RingRadii[ringnr]);      

          for ( iOme0 = iOmeMin ; iOme0 <= iOmeMax ; iOme0++) {
            iOme = iOme0 % n_ome_bins;
            if ( iOme < 0 ) iOme = iOme + n_ome_bins;
            
            // current no and max of spots in the bin  
            int iSpot = ndata[iRing][iEta][iOme];
            int maxnspot = maxndata[iRing][iEta][iOme]; 
            
            // check if the bin is already full
            if ( iSpot >= maxnspot ) { 
               maxnspot = maxnspot + 2;  //  increase size
               oldarray = data[iRing][iEta][iOme]; 
               newarray = realloc(oldarray, maxnspot * sizeof(*newarray) );
               
               if ( newarray == NULL ) {
                  printf("Memory error: memory full?\n");
                  return 1;
               }
               
               data[iRing][iEta][iOme] = newarray; 
               maxndata[iRing][iEta][iOme] = maxnspot;
            }
      
            data[iRing][iEta][iOme][iSpot] = rowno;
            (ndata[iRing][iEta][iOme])++;
         }
      }
   }

   return 0;
}




////////////////////////////////////////////////////////////////////////////////  
void
PrintBin(int iRing, int iEta, int iOmega) {
   int i;
   int noSpots = ndata[iRing][iEta][iOmega];
   
   printf("Bin: %d %d %d has %d spots:\n", iRing, iEta, iOmega, noSpots);
   for ( i = 0 ; i < noSpots ; i++) {
       int iSpot = data[iRing][iEta][iOmega][i];
       printf("  %d %d\n", i+1, iSpot);
   }
   printf("\n");

}


////////////////////////////////////////////////////////////////////////////////
// mkdir function is slightly different for linux and windows
int 
makeDir(char * fname)
{
#ifdef WIN32
	return(mkdir(fname));
#else    // for linux   and mac?                  cd                            
	return(mkdir(fname, S_IRWXU | S_IROTH | S_IXOTH));     
#endif
}



////////////////////////////////////////////////////////////////////////////////  
int
main(int argc, char *argv[]){
   printf("\n\n\t\tIndexer v4.0\nContact hsharma@aps.anl.gov in case of questions about the SHORT3DXRD project.\n\n");
   
   time_t end, start0;
   double diftotal;   
   int returncode;
   struct TParams Params;   
   int SpotIDs; // [MAX_N_SPOTS];
   int nSpotIDs;       
   char *ParamFN;   
   char fn[1024];
   
   // get command line params
   if (argc != 3) {
      printf("Supply a parameter file as argument: ie %s param.txt\n\n", argv[0]);
      exit(EXIT_FAILURE);
   }
   
   // read parameter file
   ParamFN = argv[1];
   SpotIDs = atoi(argv[2]);
   printf("Reading parameters from file: %s.\n", ParamFN);   
   returncode = ReadParams(ParamFN, &Params);
   if ( returncode != 0 ) {
      printf("Error reading params file %s\n", ParamFN );    
      exit(EXIT_FAILURE);
   }
   printf("Finished reading parameters.\n");
   
   // create outputfolder
   struct stat st;
   if (stat(Params.OutputFolder, &st) == 0) {
      printf("\nUsing existing folder for output: %s\n", Params.OutputFolder);
   }
   else {    
      returncode = makeDir(Params.OutputFolder); // for linux:, S_IRWXU | S_IROTH | S_IXOTH);
      if ( returncode != 0 ) {
         printf("Error creating folder: %s\n", Params.OutputFolder );    
         exit(EXIT_FAILURE);
      }
      printf("\nCreated folder for output: %s\n", Params.OutputFolder);      
    }   
   
   // write the params TODO windows /linux path /
   //MakeFullFileName(fn, Params.OutputFolder, "params_read.txt");
   //WriteParams(fn, Params);
   
   // generate hkls   
   printf("\nGenerating hkl's for %d rings of CellStruct = %d ...\n", Params.NrOfRings, Params.CellStruct);
   switch  ( Params.CellStruct )
   {
      case 1 : GenerateHKLsBCC(Params.RingNumbers, Params.NrOfRings, hkls, &n_hkls);
               break;
      case 2 : GenerateHKLsFCC(Params.RingNumbers, Params.NrOfRings, hkls, &n_hkls);
               break;   
      default : printf("Invalid input: invalid number for CellStruct given!\n");
             exit(EXIT_FAILURE);   
   }
   printf("No of hkl's: %d\n", n_hkls);

   // read spots
   printf("\nReading spot file: %s ...\n", Params.SpotsFileName);   
   returncode = ReadMatrixCor(Params.SpotsFileName, ObsSpotsLab, &n_spots );
   if ( returncode != 0 ) {
      printf("Error reading spot file %s\n", Params.SpotsFileName );    
      exit(EXIT_FAILURE);
   }
   printf("Finished reading spot file.\nNo of spots: %d\n\n", n_spots);

  // read ids to use for indexing
  //SpotIDs = malloc(MAX_N_SPOTS * sizeof(*SpotIDs));
  //printf("Reading IDs file: %s\n", Params.IDsFileName);
  //returncode = ReadArrayI(Params.IDsFileName, SpotIDs, &nSpotIDs);
  //if ( returncode != 0 ) {
  //printf("Error reading IDs file: %s\n", Params.IDsFileName );
  //exit(EXIT_FAILURE);
  //}
   nSpotIDs = 1;
   printf("Finished reading IDs file.\nNumber of spots to do indexing on: %d\n\n", nSpotIDs);   
   
   // add column to spot file: for each spot the distance to ideal ring
   CalcDistanceIdealRing(ObsSpotsLab, n_spots, Params.RingRadii );
         
   // put spots into bins: eta, omega, ringnr (NB ringnr 1 is stored at index 0 here!)
   printf("Binning data...\n");         
   MakeBins(ObsSpotsLab, n_spots, Params.MarginOme, Params.MarginEta, Params.StepsizeOrient, Params.RingRadii, //Params.NrOfRings,
      Params.EtaBinSize, Params.OmeBinSize);
   printf("No of bins for rings : %d\n", n_ring_bins);   
   printf("No of bins for eta   : %d\n", n_eta_bins);   
   printf("No of bins for omega : %d\n", n_ome_bins);   
   printf("Total no of bins     : %d\n\n", n_ring_bins * n_eta_bins * n_ome_bins);   
   BinStats(n_ring_bins, n_eta_bins, n_ome_bins);
   printf("Finished binning.\n\n");

   time(&start0);
   DoIndexing(SpotIDs,  nSpotIDs, Params );
   time (&end);
      
   diftotal = difftime(end, start0); 
   printf("\nTotal time elapsed [s] [min]: %.0f %.0f\n", diftotal, diftotal/60);
   
   //free(SpotIDs);
   
   return(0);
}


