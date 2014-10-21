#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iterator>
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <gsl/gsl_linalg.h>
#include <map>

using namespace std;

typedef struct clefts cleft;
struct clefts{
  float min_x,max_x,min_y,max_y,min_z,max_z;
  int width,height,depth;
  int vrtx_count;
};

struct nodeI{
  int id;
  int neibrs;
};

struct pwRun{
  string mif1;
  string mif2;
  string rnc1;
  string rnc2;
};

struct vertex {
  float coor[3];
  float ncoor[3];
  int grid[4];
  int cg[4];
  int pb[6];
  int m[6];
  int bu;
  int id;
  vector<int> env;
};

//Atom struc
typedef struct atoms atom;
struct atoms{
  float coor[3];
  float ncoor[3];
  string resn;
  string atomn;
  string chain;
  int resnb; //residue ID in PDB
  int atomnb; //atomnb in PDB
  int id;
  int bs;
  int mif;
};

//Atom struc
struct pseudoC{
  float coor[3];
  string type;
  int id;
};

typedef struct nodes node;
typedef node *pNode;
struct nodes{
  vertex* a;
  vertex* b;
  atom* ca;
  atom* cb;
  pseudoC* pa;
  pseudoC* pb;
  // float** nrgs;
  // int* pbPass;
  // int neibrs;
};

typedef struct CliqueStruct Clique;
struct CliqueStruct{
  int cg;
  vector<node> nodes;
  vector<vertex> va;
  vector<vertex> vb;
  int nbNodes;
  float tani;
  float taniX;
  float rmsd;
  float ligRMSD;
  float envD;
  float buD;
  gsl_matrix *mat_r;
  float cen_a[3];
  float cen_b[3];
  double det;
};

int *compsub = NULL;
int c,numNodes;
int numEdges=0;
int Clique_threshold=3;
int stopBk = 0;
float topT=-1.0;
float topN=-1;
int bkAll = 0;
int nCliques=0;
int maxCliques=200;
vector<Clique> cliques;

int emptOut=0;
char outH[4000];
int printDetails=0;
string pairwiseF;
string nrg_file1;
string nrg_file2;
char cmdLine[550];
char outbase[200];
char probesListFile[150];
int* probesList=NULL;
int probesListNb;
int commonInt=1;

int jttt=5;
int jtt[20][20];
int nb_of_probes=6;
int cg2=1;
char out_file[150]; //Output filename, constructed using outbase in get_info()
char outPDB[150];
char outPDB2[150];
int maxNodes=100000;
float dDist=3.0;
float ca_dDist=3.0;
float ps_dDist=1.0;
float neibr_dDist=3.0;
int cg_start=-1;
vector<vertex> mif1;
vector<vertex> mif2;
vector<atom> prot1;
vector<atom> prot2;
vector<pseudoC> pseudoL1;
vector<pseudoC> pseudoL2;
int caSize1=0;
int caSize2=0;
vector<int> list1;
vector<int> list2;
vector<int> steps;
int wrfn=0; //write sim in filename
int ss1[4];
int ss2[4];
int totalVrtx=0;
map<int,int> topCliques;
vector<atom> lig1;
vector<atom> lig2;
string rnc1;
string rnc2;
vector<pwRun> pw;

void getPairwise();
void sortArray(int *&, int nn, bool*&);
bool myfunction (nodeI,nodeI);
void AddNewClique(int, int*, int, vector<node>&);
int storeProbesList();
// bool compareSim(node ,node);
void bk(int, vector<node>&, bool*&);
void adjMat(vector<node>&, bool*&, int);
void Extend(int* old,int ne,int ce, int cg, vector<node>&, bool*&, int);
int open_file_ptr(FILE**, char*, int);
float dist3d(float[], float[]);
void clearStep(int);
void printNodes();
int read_commandline(int, char*[]);
int createVrtxVec(string, vector<vertex>&, vector<atom>&, int*, int&, vector<pseudoC>&, string, vector<atom>&);
void createNodes(int, vector<node>&, int);
void rem_spaces(char*);
int get_info(string, string);
double calcRot(vector<float>, vector<float>, float*, float*, gsl_matrix *);
double SupSVD(gsl_matrix *);
double gsl_matrix_Det3D(gsl_matrix *);
long long ConnID(int, int);
void setJTT(int);
int nam2n(string);
int samePseudo(pseudoC&,pseudoC&);