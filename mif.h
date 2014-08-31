#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <math.h> 
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iterator>

using namespace std;

#define NEW(p,type)     if ((p=(type *) malloc (sizeof(type))) == NULL) { \
        printf ("Could not allocate ptr!\n");                                    \
        exit(0);                                                        \
}

#define PI 3.14159265

//Atom struc
typedef struct atoms atom;
struct atoms{
  float x,y,z,xr,yr,zr;
  int atomId; //unique atom Id it represents
  string resn;
  string atomn;
  string chain;
  int resnb; //residue ID in PDB
  int atomnb; //atomnb in PDB
  int dir;
  int rDir;
  int h; //if its a hydrogen 0 or 1
  int bs; //if its a binding site atom
  int mif; //if it should be used to calculated MIFS
};

//Grid point struc
typedef struct vertexes vertex;
struct vertexes{
  float x,y,z;
  // float** nrg;
  // int** pb_bool;
  int* ints;
  int grid[4];
  map<string,float> env;
  int id;
};

struct atomSingleStruc{
  float* charge;
  string* res;
  string* atom; 
  int* atId;
  int count;
};

struct atPairStruc{
  float* sigma;
  float* epsilon;
  float * epsilon_FAINT;
  float* alpha_FAINT;
  float* epsilon_Custom;
  float* alpha_Custom;
};

struct atSingleStruc{
  string* name;
  int* HbD_FAINT;
  int* HbA_FAINT;
  int* aromatic_FAINT;
  int* charged_FAINT;
  int* hydrophobic_FAINT;
};

atomSingleStruc* atomSingle;
atSingleStruc* atSingle;
atPairStruc* atPair;

string cmdLine="";
string cleftFile="";
string proteinFile="";
string outBase="";
string tag="";
string chain="";
string basePath="";
string resnumc="";
float stepsize=0.5;
float maxGridDist=4.5;
float minGridDist=1.5;
float atmPbMaxDist=8.0;
float gridLigDist=2.0;
float caT=5.0;
int printDetails=0;
int printAtoms=1;
int* probesList;
float* pbDistTmax;
float* pbDistTmin;
int nbOfAts=0;
int nbOfAtoms=0;
int nbOfProbes=0;
int ss[4];
float* epsilons;
map<string,string> atomTypes;
map<string,int> eps;
map<string,int> hyd;
map<string,int> arm;
map<string,int> don;
map<string,int> acc;
map<string,int> chr;
map<string,float> minD;
map<string,float> maxD;
map<string,float> nrgT;
vector<string> probes;
vector<string> aa;

class  Protein{
  public:

    Protein(string);
    ~Protein(void);
    void readPDB(string);
    void createProteinObject();
    int getRefAtom(float&, float&, float&, string, int, string, string, int, float, float, float, string);
    vector<atom> PROTEIN;
    vector<atom> PROTEINOBJ;
    vector<float> LIGAND;
  private:
};

class  Grid{
  public:

    Grid(string, Protein&);
    ~Grid(void);
  	int readGetCleft(string, vector<atom>&, vector<float>&);
  	int generateID(int, int, int, int, int);
  	void getMinMax(string);
  	float roundCoord(float, int);
  	int inGridRes(vertex&, float);
    int smooth();
    void writeMif(vector<atom>&);
    map<int,vertex> GRID;
    vector<int> vrtxIdList;

  private:
    float min_x, min_y, min_z, max_x, max_y, max_z;
    int vrtx025,vrtx050,vrtx100,vrtx150,vrtx200;
    int width, height;
};

int readCmdLine(int, char**);
float dist_3d(float, float, float, float, float, float);
void getMif(map<int,vertex>&, vector<atom>&,vector<int>&);
int is_coord_in_cube(float,float,float,float,float,float,float);
void stripSpace(string &);
double calcNrg(vertex&, atom&, int, int&);
void getAtomRef();
void getEpsilons();
void getAtomTypes();
void getProbes();
void getaa();
