// IsoMIF is a program to identify molecular interaction field similarities between proteins
// Copyright (C) 2015 - Matthieu Chartier (under the supervision or Rafael Najmanovich)

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "mif.h"

//-----------------------------START OF MAIN-----------------------------
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/

int main(int argc, char **argv)
{
  size_t end;
  int i;

  //Reading command line arguments to store parameters, input files and output files
  if(readCmdLine(argc, argv)==24){ return(0); }

  //On vérifie s'il faut calculer les Mifs pour plusieurs protéines (-pp mon_fichier_texte)
  if(pairwiseF.compare("")!=0){ //Si un fichier contenant des paires protéine - cavité (getcleft) est défini
    getPairwise(); //on store les paires dans le vecteur px
  }else{ //S'il n'y a pas de fichier défini avec -pp
    pwRun npw;
    npw.pdbF= proteinFile;
    npw.cleftF= cleftFile;
    npw.rnc= resnumc;
    npw.ligF= ligFile;
    pw.push_back(npw);
  }

  //Reading reference files
  getAtomRef();
  getPseudoC();
  getEpsilons();
  getAtomTypes();
  getProbes();
  getaa();
  readLigFile();

  //Do all the Mif calculations run for all pairs (one or more depending on -pp)
  //On run refers to one protein PDB file and on cavity file, but we could adapt it to surface mode
  for(int pwr=0; pwr<pw.size(); pwr++){

    //Set initial values of this pair run
    proteinFile=pw[pwr].pdbF;
    cleftFile=pw[pwr].cleftF;
    resnumc=pw[pwr].rnc;
    ligFile=pw[pwr].ligF;

    //Get prefix if tag is empty (tag is the -t argument that defines an output filename base)
    if(tag.compare("")==0){
      end=proteinFile.find(".pdb");
      for(i=end-1; i>=0; i--){
        if(proteinFile.at(i)=='/'){
          tag=proteinFile.substr(i+1,end-i-1);
          break;
        }else if(i==0){
          tag=proteinFile.substr(i,end-i);
          break;
        }
      }
    }

    //Read protein PDB file and store atoms in PROTEIN vector and width,height,depth of grid
    Protein protein=Protein(proteinFile);

    //Build the grid vertices
    Grid grid=Grid(cleftFile,protein);

    //Calculate the MIFs at the vertices built
    getMif(grid.GRID,protein.PROTEIN,grid.vrtxIdList);

    //Create output files containing mif information (nrg of the probes, etc) in different format than the .mif file
    // getStats(grid.GRID,protein.LIGATOMS,grid.vrtxIdList);

    //A function that finds the closest distance to each of the 20 amino acids at each vertex and stores this distance
    //These like other descriptors can be used in the graph matching (Env is for environment)
    // getEnv(grid.GRID,protein.PROTEIN,grid.vrtxIdList);

    //Smooth is a function that looks if the potential of surrounding probes is passed the threshold and assigns this probe to the reference vrtx.
    //There is a smooth distance
    // if(smoothDist!=0) grid.smooth();

    //This function creates a vector of pseudocenters (defined in reference file).
    //A pseudocenter is a functional atom of the residues that could be important for certain interactions
    //The pseudocenters can be used in a first search step of the similarity algorithm
    //getPseudo(grid.GRID,protein.PROTEIN,grid.vrtxIdList);

    //Write the MIF and other information in the .mif output file
    grid.writeMif(protein.PROTEIN,protein.LIGATOMS);

    //Clear all vectors and maps to get ready for the next round
    tag="";
    for(int gi=0; gi<4; gi++){ ss[gi]=0; }
    for(int gi=0; gi<4; gi++){ ssm[gi]=0; }
    protein.PROTEIN.clear();
    protein.LIGAND.clear();
    protein.LIGATOMS.clear();
    grid.GRID.clear();
    grid.vrtxIdList.clear();
  }

  return(0);
}


//-------------------------------END OF MAIN-----------------------------
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/

int readCmdLine(int argc, char **argv){
  stringstream usage;
  int nb_arg;

  cout<<"Reading command line..."<<endl;

  for(nb_arg=0; nb_arg<argc; nb_arg++){
    cmdLine = cmdLine + argv[nb_arg];
    if(nb_arg<argc-1){
      cmdLine = cmdLine + " ";
    }
  }
    
  usage << "\n!---   IsoMIF   ---!\nWelcome.Bienvenue\n";
  usage << "-p <filename>        : \t Protein (PDB format)\n";
  usage << "-g <filename>        : \t Cleft of protein (GetCleft sphere format)\n";
  usage << "-surf                : \t surface mode to create the grid on the whole surface\n";
  usage << "-pp <filename>       : \t File with highthroughput mif calculations\n";
  usage << "-s                   : \t Which grid resolution to consider [0.25, 0.5, 1.0 or 2.0] default is 0.5 Angstrom\n";
  usage << "-d <distance in Ang.>: \t atom-probe distance threshold. (default: "<< atmPbMaxDist <<" Angstrom)\n";
  usage << "-c <chain>           : \t PDB chains to consider for grid and energy calculations\n";
  usage << "-w <distance in Ang.>: \t minimum distance of the probes to the atoms (default: "<< minGridDist <<" Angstrom)\n";
  usage << "-v <distance in Ang.>: \t maximum distance of the probes to the atoms (default: "<< maxGridDist <<" Angstrom)\n";
  usage << "-pr                  : \t print details of energy calculations for each probe-atom pair\n";
  usage << "-o                   : \t output directory\n";
  usage << "-b                   : \t burriedness level (0-14) [8]\n";
  usage << "-m                   : \t smooth\n";
  usage << "-ff                  : \t force field to use\n";
  usage << "-mat                 : \t epsilon matrix to use\n";
  usage << "-pb                  : \t probes threshold file to use\n";
  usage << "-a                   : \t angle threshold\n";
  usage << "-t                   : \t filename (tag) [pdb name by default]\n";
  usage << "-lf                  : \t ligand file with atom types\n";
  usage << "-dpv                 : \t distance probe-Vrtx to consider a vrtx to be true-positive\n";
  usage << "-sf                  : \t stats output file\n";
  usage << "-l                   : \t RESNUMC of the ligand from which to crop the grid\n";
  usage << "-r                   : \t maximum distance between the grid and the ligand or residues atoms\n";
  usage << "-x                   : \t do not write atoms in mif file\n";
  usage << "-og                  : \t dir for grid file\n";
  usage << "-h                   : \t help menu\n";

  if(argc<5){
    cout << endl << "Missing obligatory arguments:" << usage.str() << endl;
    return(24);
  }

  basePath=argv[0];
  unsigned found = basePath.find_last_of("/\\");
  basePath=basePath.substr(0,found) + '/';

  for(nb_arg=1; nb_arg<argc; nb_arg++){
    if(strcmp(argv[nb_arg],"-h")==0){
      cout<< usage.str() << endl;
      return(24);
    }
    if(strcmp(argv[nb_arg],"-p")==0){
      proteinFile=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-g")==0){
      cleftFile=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-x")==0){
      printAtoms=0;
    }
    if(strcmp(argv[nb_arg],"-o")==0){
      outBase=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-ff")==0){
      ff=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-sf")==0){
      statsF=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-c")==0){
      chain=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-l")==0){
      resnumc=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-lf")==0){
      ligFile=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-mat")==0){
      matrixF=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-pb")==0){
      probesF=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-m")==0){
      sscanf(argv[nb_arg+1], "%d", &smoothDist);
    }
    if(strcmp(argv[nb_arg],"-v")==0){
      sscanf(argv[nb_arg+1], "%f", &maxGridDist);
    }
    if(strcmp(argv[nb_arg],"-w")==0){
      sscanf(argv[nb_arg+1], "%f", &minGridDist);
    }
    if(strcmp(argv[nb_arg],"-a")==0){
      sscanf(argv[nb_arg+1], "%f", &angThresh);
    }
    if(strcmp(argv[nb_arg],"-r")==0){
      sscanf(argv[nb_arg+1], "%f", &gridLigDist);
    }
    if(strcmp(argv[nb_arg],"-d")==0){
      sscanf(argv[nb_arg+1], "%f", &atmPbMaxDist);
    }
    if(strcmp(argv[nb_arg],"-dpv")==0){
      sscanf(argv[nb_arg+1], "%f", &distpbV);
    }
    if(strcmp(argv[nb_arg],"-s")==0){
      sscanf(argv[nb_arg+1], "%f", &stepsize);

      if(stepsize == 0.5){
        zip=3;
      }else if(stepsize == 1.0){
        zip=2;
      }else if(stepsize == 1.5){
        zip=1;
      }else if(stepsize == 2.0){
        zip=0;
      }
    }
    if(strcmp(argv[nb_arg],"-b")==0){
      sscanf(argv[nb_arg+1], "%d", &bul);
    }

    if(strcmp(argv[nb_arg],"-og")==0){
      outGridBase=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-pr")==0){
      printDetails=1;
    }
    if(strcmp(argv[nb_arg],"-t")==0){
      tag=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-pp")==0){
      pairwiseF=argv[nb_arg+1];
    }
    if(strcmp(argv[nb_arg],"-surf")==0){
      surf =1;
    }
  }

  if(chain.compare("")==0){
    chain="none";
  }
  if(outBase.compare("")==0){
    outBase="./";
  }else{
    outBase = outBase + "/";
  }

  //Remove the alternate location
  resnumcShort=resnumc.substr(0,resnumc.length()-1);

  cout<<endl;
  cout<< "Protein file: "<< proteinFile <<endl;
  cout<< "GetCleft file: "<< cleftFile <<endl;
  cout<< "Output file base: "<< outBase <<endl;
  cout<< "Chain: "<< chain <<endl;
  cout<< "Tag: " << tag << endl;
  cout<< "Burriedness level: " << bul << endl;
  cout<< "Stepsize: "<<stepsize<<endl;
  cout<< "Gridstep: "<<gridStep<<endl;
  cout<< "MinGridDist: "<< minGridDist <<endl;
  cout<< "MaxGridDist: "<< maxGridDist <<endl;
  cout<< "AtmPbMaxDist: "<< atmPbMaxDist <<endl;
  cout<< "RESNUMC: "<< resnumc <<endl;
  cout<< "RESNUMC short: "<< resnumcShort <<endl;
  cout<< "gridLigDist: "<< gridLigDist <<endl;

  //Square some distances so we don't need to do the square root in distance calculations
  gridLigDist=gridLigDist*gridLigDist;
  maxGridDist=maxGridDist*maxGridDist;
  minGridDist=minGridDist*minGridDist;

  return(0);
}

Protein::~Protein(void){

}

Protein::Protein(string filename){
  min_x=1000.00;
  min_y=1000.00;
  min_z=1000.00;
  max_x=-1000.00;
  max_y=-1000.00;
  max_z=-1000.00;

  //Read PDB file of the protein and store atoms
  readPDB(filename);

  //Get reference atom for each atom that will be used to calculate angles
  getAtomDir();
}

void Protein::readPDB(string filename){
  ifstream ifs;
  ofstream ofs;
  string line;
  string fields[13];
  string thisresnumc;
  string thisresnumcShort;
  float x,y,z;
  int atomnb;
  int resnb;
  atom* atm=NULL;
  string outFileName=outBase + tag + "_cpy.pdb";
  size_t found;
  float minx,miny,minz,maxx,maxy,maxz;

  cout<<endl<<"Reading PDB file..."<<endl;

  //Opening PDB output file that will contain a copy of the PDB structure
  //This file is always in the same directory as the .mif file and is used in mifView.pl to create a object of the structure in PyMOL
  ofs.open(outFileName.c_str());

  ifs.open(filename.c_str());
  if(!ifs.is_open()){ cout << "could not read "<< filename << endl; }
  while(ifs.good()){

    getline(ifs,line);

    if(line.compare(0,3,"END") == 0){ break; }

    //Print copy of the PDB
    ofs << line << endl; 
        
    if(line.compare(0,6,"ATOM  ") == 0 || line.compare(0,6,"HETATM") == 0){
      fields[0] = line.substr(0,6);   // ATOM/HETATM field
      fields[1] = line.substr(6,5);   // atom number
      fields[2] = line.substr(12,4);  // atom name
      fields[3] = line.substr(16,1);  // alternative
      fields[4] = line.substr(17,3);  // residue name
      fields[5] = line.substr(21,1);  // residue chain
      fields[6] = line.substr(22,4);  // residue number
      fields[8] = line.substr(30,8);  // x-coord
      fields[9] = line.substr(38,8);  // y-coord
      fields[10] = line.substr(46,8); // z-coord
      
      //Remove spaces
      stripSpace(fields[1]);
      stripSpace(fields[2]);
      stripSpace(fields[3]);
      stripSpace(fields[4]);
      stripSpace(fields[5]);
      stripSpace(fields[6]);

      atomnb=atoi(fields[1].c_str());
      resnb=atoi(fields[6].c_str());
      x=atof(fields[8].c_str());
      y=atof(fields[9].c_str());
      z=atof(fields[10].c_str());

      if(fields[4].compare("HOH")==0){ continue; } //Continue if its a water molecule line

      atom atm;
      atm.x=x;
      atm.y=y;
      atm.z=z;
      atm.chain=fields[5];
      atm.alt=fields[3];
      atm.atomnb=atomnb;
      atm.resnb=resnb;
      atm.atomn=fields[2];
      atm.resn=fields[4];
      atm.bs=0;
      atm.mif=1;
      atm.h=0;

      if(atm.alt.compare("")==0){ atm.alt="-"; } //If there is no alterante location for this atom set alt to '-'
      if(atm.chain.compare("")==0){ atm.chain="-"; } //If there is no chain for this atom set chain to '-'

      //See if the atom in this line is a hydrogen (atm.h = 1)
      try{
        fields[11] = line.substr(76,2); //Get atom symbol of this atom
        stripSpace(fields[11]);
        found = fields[11].find("H"); //If this atom symbol contains an H
        if(found!=string::npos) atm.h=1; //Tag this atom as a hydrogen atom
      }catch(...){
        // cout<<"No atom symbol"<<endl;
      }
      
      //Check conditions and set atm.mif to 0 if necessary
      //atm.mif == 0 when the atom must not be considered to build the grid or calculate the Mif
      if(atomTypes.find(atm.resn+"_"+atm.atomn) == atomTypes.end()){ atm.mif=0; } //If theres no atom type in our files
      if(fields[5].compare(chain)!=0 && chain.compare("none")!=0){ atm.mif=0; } //If its not in the PDB chain we want
      if(line.compare(0,6,"HETATM") == 0){ atm.mif=0; } //If its an HETATM atom (bound ligand, ions, water molecules)
      if((fields[3].compare("A")!=0) && (fields[3].compare("")!=0)){ atm.mif=0; } //If it has an alternate location and its not A

      //If we defined an anchor amino acid or residue (-l resnumc) used to crop the grid around it
      //we must store the anchor atoms in LIGAND
      if(resnumc.compare("")!=0){

        //Get long and short resnumc (with and without alt location)
        stringstream sss;
        sss << atm.resnb; //We put the number in a stringstream to convert it to a string easily
        thisresnumc = atm.resn + sss.str() + atm.chain + atm.alt; //long resnumc
        thisresnumcShort = atm.resn + sss.str() + atm.chain; //short resnumc
        stripSpace(thisresnumc);

        //If its not a hydrogen atom and the resnumc defined in the cmd line (with -l argument) is the same as this atom, push as LIGAND atom
        if(atm.h==0 && (resnumc.compare(thisresnumc)==0 || (resnumcShort.compare(thisresnumcShort)==0 && atm.alt.compare("-")==0))){
          LIGAND.push_back(atof((line.substr(30,8).c_str())));
          LIGAND.push_back(atof((line.substr(38,8).c_str())));
          LIGAND.push_back(atof((line.substr(46,8).c_str())));
          atm.at=ligAt[atm.atomn];
          LIGATOMS.push_back(atm);
        }
      }

      //Round min/max X,Y and Z coords
      minx=roundCoord(x-0.5,0);
      miny=roundCoord(y-0.5,0);
      minz=roundCoord(z-0.5,0);
      maxx=roundCoord(x+0.5,1);
      maxy=roundCoord(y+0.5,1);
      maxz=roundCoord(z+0.5,1);

      //Update min/max X,Y and Z intervals
      if(minx<min_x){ min_x=minx; }
      if(miny<min_y){ min_y=miny; }
      if(minz<min_z){ min_z=minz; }
      if(maxx>max_x){ max_x=maxx; }
      if(maxy>max_y){ max_y=maxy; }
      if(maxz>max_z){ max_z=maxz; }

      //Push the atom in the protein vector
      PROTEIN.push_back(atm);
    }
  }
  ifs.close();
  ofs.close();

  //Add buffer layer (maximum distance allowed between an atom and a grid point + stepsize) to the min/max X,Y and Z
  min_x-=maxGridDist+gridStep; min_y-=maxGridDist+gridStep; min_z-=maxGridDist+gridStep; max_x+=maxGridDist+gridStep; max_y+=maxGridDist+gridStep; max_z+=maxGridDist+gridStep;
  
  //Calculate width, height and depth of the mother grid
  width=(int)(((max_x-min_x)/gridStep)+1.0);
  height=(int)(((max_y-min_y)/gridStep)+1.0);
  depth=(int)(((max_z-min_z)/gridStep)+1.0);

  //Print values to the terminal
  cout<<"PROTEIN min/max values: minx: "<< min_x<< " miny: "<< min_y<< " minz: "<<min_z<<" maxx: "<<max_x<<" maxy: "<<max_y<<" maxz: "<<max_z<<endl;
  cout<<"PROTEIN Width "<<width<<" Height "<< height << endl;
  cout<<"Number of anchor atoms: "<<LIGAND.size()<<endl;

}

void Protein::getAtomDir(){
  int i,found,ring,needRef,rDir;
  float xr,yr,zr;
  string tatomn;
  string tatomn2;

  for(i=0; i<PROTEIN.size(); i++){ //For each atom of the protein
    if(PROTEIN[i].mif==1){ //If its an atom we want to consider for the Mif calculations
      ring=0;
      found=0;
      needRef=0;
      if(atomRef[PROTEIN[i].resn+"_"+PROTEIN[i].atomn].r1.compare("")!=0){ // If there is a reference atom based on our input files
        needRef=1;
        tatomn=atomRef[PROTEIN[i].resn+"_"+PROTEIN[i].atomn].r1;
        rDir=atomRef[PROTEIN[i].resn+"_"+PROTEIN[i].atomn].rDir;
        if(atomRef[PROTEIN[i].resn+"_"+PROTEIN[i].atomn].r2.compare("")!=0){ // If its a ring atom
          tatomn2=atomRef[PROTEIN[i].resn+"_"+PROTEIN[i].atomn].r2;
          ring=atomRef[PROTEIN[i].resn+"_"+PROTEIN[i].atomn].ring;
        }
      }
      if(needRef==1){
        if(getRefAtom(xr, yr, zr, PROTEIN[i].resn, PROTEIN[i].resnb, tatomn, tatomn2, ring, PROTEIN[i].x, PROTEIN[i].y, PROTEIN[i].z, PROTEIN[i].chain)==1){
          PROTEIN[i].xr=xr;
          PROTEIN[i].yr=yr;
          PROTEIN[i].zr=zr;
          PROTEIN[i].dir=1;
          PROTEIN[i].rDir=rDir;
        }
      }else{
        PROTEIN[i].xr=0.0;
        PROTEIN[i].yr=0.0;
        PROTEIN[i].zr=0.0;
        PROTEIN[i].dir=0;
      }
    }
  }
  cout<<"Protein has "<< PROTEIN.size()<<" atoms"<<endl<<endl;
}

int Protein::getRefAtom(float& xr, float& yr, float& zr, string tresn, int tresnb, string tatomn, string tatomn2, int ring, float x, float y, float z, string chain){
  int i;
  int found=0;
  int foundr1=0;
  int foundr2=0;
  float xr1,yr1,zr1,xr2,yr2,zr2;
  float ax,ay,az,bx,by,bz;
  float dist,rDist,rpDist,angle;

  if(ring==0){ //This atom is not part of a ring
    for(i=0; i<PROTEIN.size(); i++){
      if(PROTEIN[i].resn.compare(tresn)==0 && PROTEIN[i].atomn.compare(tatomn)==0 && PROTEIN[i].resnb==tresnb && PROTEIN[i].chain.compare(chain)==0){
        found=1;
        xr=PROTEIN[i].x;
        yr=PROTEIN[i].y;
        zr=PROTEIN[i].z;
      }
      if(found==1){ break; }
    }
  }else{ //If its part of a ring
    for(i=0; i<PROTEIN.size(); i++){

      //Loop to find reference atom 1
      if(PROTEIN[i].resn.compare(tresn)==0 && PROTEIN[i].atomn.compare(tatomn)==0 && PROTEIN[i].resnb==tresnb && PROTEIN[i].chain.compare(chain)==0){
        foundr1=1;
        xr1=PROTEIN[i].x;
        yr1=PROTEIN[i].y;
        zr1=PROTEIN[i].z;
      }
      //Loop to find reference atom 2
      if(PROTEIN[i].resn.compare(tresn)==0 && PROTEIN[i].atomn.compare(tatomn2)==0 && PROTEIN[i].resnb==tresnb && PROTEIN[i].chain.compare(chain)==0){
        foundr2=1;
        xr2=PROTEIN[i].x;
        yr2=PROTEIN[i].y;
        zr2=PROTEIN[i].z;
      }

      //If we found both reference atoms get out of the for loop
      if(foundr1==1 && foundr2==1){
        break;
      }
    }

    //If we found both reference atom
    if(foundr1==1 && foundr2==1){
      found=1;
      if(ring==1){ //Ring type 1 means we must find a reference atom to calculate angle with a vector equidistant from adjacent atoms (nitrogen of histidine)
        ax=x-xr1;
        ay=y-yr1;
        az=z-zr1;
        bx=x-xr2;
        by=y-yr2;
        bz=z-zr2;
        xr=((ay*bz)-(by*az))+x;
        yr=((bx*az)-(ax*bz))+y;
        zr=((ax*by)-(ay*bx))+z;
        // dist=dist_3d(x,y,z,xr1,yr1,zr1);
        // rDist=dist_3d(x,y,z,xr,yr,zr);
        // rpDist=dist_3d(xr1,yr1,zr1,xr,yr,zr);
        // angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
        // angle=acos(angle)* 180.0 / PI;
      }else if(ring==2){ //Ring type 2 means we must find a reference atom to calculate angle between probe and plane of the ring
        xr=(xr1-x+xr2-x)+x;
        yr=(yr1-y+yr2-y)+y;
        zr=(zr1-z+zr2-z)+z;
        // dist=dist_3d(x,y,z,xr1,yr1,zr1);
        // rDist=dist_3d(x,y,z,xr,yr,zr);
        // rpDist=dist_3d(xr1,yr1,zr1,xr,yr,zr);
        // angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
        // angle=acos(angle)* 180.0 / PI;
        // cout<<"angle 1 "<<angle<<endl;
        // dist=dist_3d(x,y,z,xr2,yr2,zr2);
        // rDist=dist_3d(x,y,z,xr,yr,zr);
        // rpDist=dist_3d(xr2,yr2,zr2,xr,yr,zr);
        // angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
        // angle=acos(angle)* 180.0 / PI;
        // cout<<"angle 2 "<<angle<<endl;
        // cout << x<<" "<<y<<" "<<z<<" ref1 "<<xr1<<" "<<yr1<<" "<<zr1<<" ref2 "<<xr2<<" "<<yr2<<" "<<zr2<< " vec "<<xr<<" "<<yr<<" "<<zr<<endl;
      }
    }
  }

  return(found);
}

Grid::Grid(string filename, Protein& prot){
  vrtx050=0;
  vrtx100=0;
  vrtx150=0;
  vrtx200=0;

  if(filename.compare("")==0){
    //Build the grid everywhere around the protein and classify them as protein vertices (p=1) or probe vertices (p=0)
    buildGrid(prot);

    //Calculate burriedness of each vertex (value of 0 to 14), 0 being completely solvent exposed and 14 completely burried
    getBuriedness();
  }else{

    //Create grid using the getcleft spheres
    readGetCleft(filename, prot.PROTEIN, prot.LIGAND);

    //If we want to filter for burriedness we must find the protein vertices (vrtx.p==1)
  }
}

Grid::~Grid(void){

}

int Grid::generateID(int w, int h, int x, int y, int z){
  int id;
  id=z * (w * h) + y * w + x;
  return(id);
}

int Grid::buildGrid(Protein& prot){ //Building the grid using the whole protein
  int i=0;
  int pggrid2=0;

  cout<< "Building Grid..."<<endl;

  //Go through our min/max X,Y and Z intervals to build the mother GRID
  map<int,vertex>::iterator it;
  int id;
  for(float x=min_x; x<=max_x; x+=stepsize){
    for(float y=min_y; y<=max_y; y+=stepsize){
      for(float z=min_z; z<=max_z; z+=stepsize){

        //Check distance to closest protein atom
        int pg=0;
        float minDist=10000.0;
        for(i=0; i<prot.PROTEIN.size(); i++){
          if(prot.PROTEIN[i].mif==0) continue; //Skip if its not an atom that we consider (respects the filters in readPDB)
          float d=((abs(x-prot.PROTEIN[i].x)*abs(x-prot.PROTEIN[i].x))+(abs(y-prot.PROTEIN[i].y)*abs(y-prot.PROTEIN[i].y))+(abs(z-prot.PROTEIN[i].z)*abs(z-prot.PROTEIN[i].z)));
          if(d<minDist){ minDist=d; }
        }

        if(minDist<minGridDist){ //Its too close to the protein atoms so its a protein vertex
          pg=1;
        }else if(minDist > maxGridDist){ //Its too far from the protein so its skipped
          continue;
        }

        //Generate grid vertex ID
        id=generateID(width,height,(int)((x-min_x)/gridStep)+1,(int)((y-min_y)/gridStep)+1,(int)((z-min_z)/gridStep)+1);
        it = GRID.find(id);

        if(it == GRID.end()){ //If we never added this point to the vector
          vertex vrtx;
          vrtx.x=x;
          vrtx.y=y;
          vrtx.z=z;

          if(pg==1){ //If its a grid vertex that represents a protein vertex
            vrtx.p=1;
            vrtx.id=uID;
            uID++;
          }else{ //If its a grid vertex where we'll calculate a MIF
            vrtx.p=0;

            //Calculate in which grid resolution this vertex is
            if(inGridRes(vrtx,2.0)==1){
                vrtx.grid[0]=1;
            }else{ vrtx.grid[0]=0; }

            if(inGridRes(vrtx,1.5)==1){
              vrtx.grid[1]=1;
            }else{ vrtx.grid[1]=0; }

            if(inGridRes(vrtx,1.0)==1){
              vrtx.grid[2]=1;
            }else{ vrtx.grid[2]=0; }

            if(inGridRes(vrtx,0.5)==1){
              vrtx.grid[3]=1;
            }else{ vrtx.grid[3]=0; }

            //Skip this vertex if its not in the grid resolution we want (defined by argument -s)
            if(vrtx.grid[zip]!=1 && zip!=-1) continue;

            //We must remove grid vertex that are too far from the anchor atoms
            if(resnumc.compare("")!=0){
              float minDist=10000.0;

              //Loop to find the distance to the closest anchor atom and store it in minDist
              for(int p=0; p < prot.LIGAND.size(); p+=3){
                float d=((abs(x-prot.LIGAND[p])*abs(x-prot.LIGAND[p]))+(abs(y-prot.LIGAND[p+1])*abs(y-prot.LIGAND[p+1]))+(abs(z-prot.LIGAND[p+2])*abs(z-prot.LIGAND[p+2])));
                if(d < minDist){minDist = d;}
              }
              //If minDist is bigger than the threshold distance (-r argument)
              if(minDist > gridLigDist) continue; //Skip to next vertex
            }

            //Initialize probe interactions, nrg potential and angle values of this vertex
            int zi=0;
            float zf=0.0;
            for(int i=0; i<nbOfProbes; i++){
              vrtx.ints.push_back(zi); vrtx.nrgs.push_back(zf); vrtx.angles.push_back(zf);
            }

            //Initialize the environment vector for each amino acid
            // for(int i=0; i<aa.size(); i++){ vrtx.env[aa[i]]=1000.0; }

            if(vrtx.grid[0]==1){ vrtx200++; }
            if(vrtx.grid[1]==1){ vrtx150++; }
            if(vrtx.grid[2]==1){ vrtx100++; }
            if(vrtx.grid[3]==1){ vrtx050++; }

            vrtx.id=uID;
            uID++;
          }

          //Insert the vertex in the grid
          GRID.insert(pair<int,vertex>(id,vrtx));
          vrtxIdList.push_back(id);
        }
      }
    }
  }

  cout<<"Grid points [2.0] "<<vrtx200<<" [1.5] "<<vrtx150<<" [1.0] "<<vrtx100<<" [0.5] "<<vrtx050<<"."<<endl;

  return(0);
}

int Grid::readGetCleft(string filename, vector<atom>& protVec, vector<float>& ligVec){
  ifstream ifs;
  string line="";
  string fields[4];
  map<int,vertex>::iterator it;
  vector<atom>::iterator pit;
  int i=0;
  int newv=0;
  int id;
  float x,y,z,rad,nx,ny,nz;
  float minx,miny,minz,maxx,maxy,maxz;
  float dist,minDist;

  ifs.open(filename.c_str());

  if(!ifs.is_open()){ 
    cout << "could not read "<< filename << endl;
  }

  cout<<"Reading cleft file..."<<endl;
  while(ifs.good()){
    getline(ifs,line);
        
    if(line.compare(0,6,"ATOM  ") == 0 || line.compare(0,6,"HETATM") == 0){
      // cout<< line<<endl;
      fields[0] = line.substr(30,8);   // x-coord
      fields[1] = line.substr(38,8);   // y-coord
      fields[2] = line.substr(46,8);   // z-coord
      fields[3] = line.substr(60,6);   // bfactor
      x=atof(fields[0].c_str());
      y=atof(fields[1].c_str());
      z=atof(fields[2].c_str());
      rad=atof(fields[3].c_str());
      rad+=1.0;

      minx=roundCoord(x-rad,0);
      miny=roundCoord(y-rad,0);
      minz=roundCoord(z-rad,0);
      maxx=roundCoord(x+rad,1);
      maxy=roundCoord(y+rad,1);
      maxz=roundCoord(z+rad,1);

      for(nx=minx; nx<=maxx; nx+=gridStep){
        for(ny=miny; ny<=maxy; ny+=gridStep){
          for(nz=minz; nz<=maxz; nz+=gridStep){

            //Generate grid vertex ID
            id=generateID(width,height,(int)((nx-min_x)/gridStep)+1,(int)((ny-min_y)/gridStep)+1,(int)((nz-min_z)/gridStep)+1);
            it = GRID.find(id);

            //Skip if already exists
            if(it != GRID.end()) continue;

            //Skip if not within the sphere
            if(((abs(nx-x)*abs(nx-x))+(abs(ny-y)*abs(ny-y))+(abs(nz-z)*abs(nz-z)))>(rad*rad)) continue;

            if(resnumc.compare("")!=0){ //If we must remove grid points too far from the anchor atoms
              minDist=10000.0;
              for(i=0; i<ligVec.size(); i+=3){
                dist=((abs(nx-ligVec.at(i))*abs(nx-ligVec.at(i)))+(abs(ny-ligVec.at(i+1))*abs(ny-ligVec.at(i+1)))+(abs(nz-ligVec.at(i+2))*abs(nz-ligVec.at(i+2))));
                // cout<< ligVec.at(i) << " "<< ligVec.at(i+1) << " "<< ligVec.at(i+2) << " grid: "<<nx<<" "<<ny<<" "<<nz<< " dist:"<< dist<< endl;
                if(dist<minDist) minDist=dist;
              }
              //Skip to next grid intersection if too far from the ligand
              if(minDist>gridLigDist){
                // cout<<nx<<" "<<ny<<" "<<nz<<"to far "<<minDist<<endl;
                continue;
              }
            }

            vertex vrtx;
            vrtx.x=nx;
            vrtx.y=ny;
            vrtx.z=nz;
            vrtx.bu=0;
            vrtx.p=0;
            vrtx.modulo=id;

            for(int i=0; i<4; i++){ vrtx.grid[i]=0; }

            int zi=0;
            float zf=0.0;
            for(int i=0; i<nbOfProbes; i++){
              vrtx.ints.push_back(zi); vrtx.nrgs.push_back(zf); vrtx.angles.push_back(zf);
            }

            //Calculate in which grid resolution this vertex is
            if(inGridRes(vrtx,2.0)==1){
                vrtx.grid[0]=1;
            }else{ vrtx.grid[0]=0; }

            if(inGridRes(vrtx,1.5)==1){
              vrtx.grid[1]=1;
            }else{ vrtx.grid[1]=0; }

            if(inGridRes(vrtx,1.0)==1){
              vrtx.grid[2]=1;
            }else{ vrtx.grid[2]=0; }

            if(inGridRes(vrtx,0.5)==1){
              vrtx.grid[3]=1;
            }else{ vrtx.grid[3]=0; }

            vrtx.id=uID;
            uID++;
            newv++;

            GRID.insert(pair<int,vertex>(id,vrtx));
          }
        }
      }
    }
  }
  ifs.close();
  cout<<"Grid points added: "<<newv<<endl;

  int tooClose=0;
  int tooFar=0;
  int badRes=0;
  //Remove points too close or too far for the protein
  cout<<endl<<"Chipping grid..."<<endl;
  int chipped=0;
  it=GRID.begin();
  while(it!=GRID.end()){

    if(it->second.p==1){ //If its a protein vertex
      ++it;
      continue;
    }

    //Get distance to the closest protein atom
    minDist=10000.0;
    for(pit=protVec.begin(); pit!=protVec.end(); ++pit){
      if((*pit).h==1) continue;
      if((*pit).mif==0) continue;
      dist=((abs(it->second.x-(*pit).x)*abs(it->second.x-(*pit).x))+(abs(it->second.y-(*pit).y)*abs(it->second.y-(*pit).y))+(abs(it->second.z-(*pit).z)*abs(it->second.z-(*pit).z)));
      if(dist<minDist) minDist=dist;
    }

    // cout<<minDist<<" "<<it->second.grid[0]<<" "<<it->second.grid[1]<<" "<<it->second.grid[2]<<" "<<it->second.grid[3]<<endl;

    if(minDist<minGridDist){
      GRID.erase(it++);
      chipped++;
      tooClose++;
    }
    else
    if(minDist > maxGridDist){
      GRID.erase(it++);
      chipped++;
      tooFar++;
    }
    else
    if(it->second.grid[zip]==0 && zip!=-1){
      GRID.erase(it++);
      chipped++;
      badRes++;
    }
    else{
      if(inGridRes(it->second,2.0)==1) vrtx200++;
      if(inGridRes(it->second,1.5)==1) vrtx150++;
      if(inGridRes(it->second,1.0)==1) vrtx100++;
      if(inGridRes(it->second,0.5)==1) vrtx050++;
      vrtxIdList.push_back(it->second.modulo);
      ++it;
    }
  }

  cout<<"Chipped "<<chipped<<" grid points."<<endl;
  cout<<"Too close "<<tooClose<<" too far "<<tooFar<<" bad res "<<badRes<<endl;
  cout<<"Grid points [2.0] "<<vrtx200<<" [1.5] "<<vrtx150<<" [1.0] "<<vrtx100<<" [0.5] "<<vrtx050<<"."<<endl;
  return(0);
}

 void Grid::getBuriedness(){
   map<int,vertex>::iterator it;
   int id;
   int nbu=0;

   //Get buriedness of each grid point
   cout<<endl<<"Getting buriedness..."<<endl;
   map<int,vertex>::iterator m=GRID.begin();
   while(m!=GRID.end()){
     if(m->second.p==1){ //If its a protein vertex, put its burriedness to 0, meaning completely burried.
       m->second.bu=0;
       ++m;
       continue;
     }
     int bu=0;

     //position of vertex in modulo
     int xi=(int)((m->second.x-min_x)/gridStep)+1;
     int yi=(int)((m->second.y-min_y)/gridStep)+1;
     int zi=(int)((m->second.z-min_z)/gridStep)+1;

     //modulo limits in the 6 directions: xl (x-left), xr (x-right), yd (y-down), yu (y-up), zb (z-back), zf (z-front)
     int xl=(int)((m->second.x-min_x)/gridStep);
     int xr=(int)((max_x-m->second.x)/gridStep);
     int yd=(int)((m->second.y-min_y)/gridStep);
     int yu=(int)((max_y-m->second.y)/gridStep);
     int zb=(int)((m->second.z-min_z)/gridStep);
     int zf=(int)((max_z-m->second.z)/gridStep);

     // cout<<xi<<" "<<yi<<" "<<zi<<endl;
     // cout<<min_x<<" | "<<xl<<" "<<m.x<<" "<<xr<<" | "<<max_x<<" || "<<width<<endl;
     // cout<<min_y<<" | "<<yd<<" "<<m.y<<" "<<yu<<" | "<<max_y<<" || "<<height<<endl;
     // cout<<min_z<<" | "<<zb<<" "<<m.z<<" "<<zf<<" | "<<max_z<<" || "<<depth<<endl;

     //Going through the 6 directions to see if a protein vertex is in the way
     // cout<<endl<<"going xl"<<endl;
     int move=0;
     for(int tx=xi-1; tx>=1; tx--){
       move++;
       if(move>buD) break;
       id=generateID(width,height,tx,yi,zi);
       it = GRID.find(id);
       if(it != GRID.end()){
         if(GRID[id].p==1){
           bu++;
           break;
         }
       }
     }

     // cout<<endl<<"going xr"<<endl;
     move=0;
     for(int tx=xi+1; tx<=width; tx++){
       move++;
       if(move>buD) break;
       id=generateID(width,height,tx,yi,zi);
       it = GRID.find(id);
       if(it != GRID.end()){
         if(GRID[id].p==1){
           bu++;
           break;
         }
       }
     }

     // cout<<endl<<"going yd"<<endl;
     move=0;
     for(int ty=yi-1; ty>=1; ty--){
       move++;
       if(move>buD) break;
       id=generateID(width,height,xi,ty,zi);
       // cout<<xi<<" "<<ty<<" "<<zi<<endl;
       it = GRID.find(id);
       if(it != GRID.end()){
         if(GRID[id].p==1){
           bu++;
           break;
         }
       }
     }
     // cout<<endl<<"going yu"<<endl;
     move=0;
     for(int ty=yi+1; ty<=height; ty++){
       move++;
       if(move>buD) break;
       id=generateID(width,height,xi,ty,zi);
       // cout<<xi<<" "<<ty<<" "<<zi<<endl;
       it = GRID.find(id);
       if(it != GRID.end()){
         if(GRID[id].p==1){
           bu++;
           break;
         }
       }
     }
     // cout<<endl<<"going zb"<<endl;
     move=0;
     for(int tz=zi-1; tz>=1; tz--){
       move++;
       if(move>buD) break;
       // cout<<xi<<" "<<yi<<" "<<tz<<endl;
       id=generateID(width,height,xi,yi,tz);
       it = GRID.find(id);
       if(it != GRID.end()){
         if(GRID[id].p==1){
           bu++;
           break;
         }
       }
     }
     // cout<<endl<<"going zf"<<endl;
     move=0;
     for(int tz=zi+1; tz<=depth; tz++){
       move++;
       if(move>buD) break;
       // cout<<xi<<" "<<yi<<" "<<tz<<endl;
       id=generateID(width,height,xi,yi,tz);
       it = GRID.find(id);
       if(it != GRID.end()){
         if(GRID[id].p==1){
           bu++;
           break;
         }
       }
     }

     //Going through the 8 diagonals to see if a protein vertex is in the way
     int tx=xi;
     int ty=yi;
     int tz=zi;
     int flag=0;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx--; ty++; tz++;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }

     flag=0;
     tx=xi; ty=yi; tz=zi;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx--; ty++; tz--;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }

     flag=0;
     tx=xi; ty=yi; tz=zi;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx++; ty++; tz--;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }

     flag=0;
     tx=xi; ty=yi; tz=zi;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx++; ty++; tz++;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }

     flag=0;
     tx=xi; ty=yi; tz=zi;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx--; ty--; tz++;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }

     flag=0;
     tx=xi; ty=yi; tz=zi;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx--; ty--; tz--;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }

     flag=0;
     tx=xi; ty=yi; tz=zi;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx++; ty--; tz--;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }

     flag=0;
     tx=xi; ty=yi; tz=zi;
     move=0;
     while(flag==0){
       move++;
       if(move>buD) break;
       tx++; ty--; tz++;
       int diag=getDiag(tx,ty,tz,flag);
       if(diag){
         bu++;
       }
     }
      if(bu>bul){
        nbu++;
      }

     if(bu<=bul){
       m->second.bu=bu;
       id=generateID(width,height,xi,yi,zi);
       vrtxIdList.push_back(id);
       ++m;
     }else{
       GRID.erase(m++);
       nbu++;
     }
   }
   cout<<"Removed "<<nbu<<" vertex with a burriedness threshold <=" << bul <<endl;
 }

int Grid::getDiag(int tx, int ty, int tz, int& flag){
  int id=generateID(width,height,tx,ty,tz);

  map<int,vertex>::iterator it = GRID.find(id);
  if(it != GRID.end()){
    if(GRID[id].p==1){
      flag=1;
      return(1);
    }else{
      return(0);
    }
  }else{
    flag=1;
    return(0);
  }
}

float roundCoord(float number, int min_or_max){ 
  //min_or_max (1=max, 0=min)
  float quotient;
  int rounded;
  float new_coord;
  quotient=number/gridStep;
  rounded=(int)quotient;
  new_coord=rounded*gridStep;

  if(min_or_max==1){ // max value
    if(new_coord<number || fabs(new_coord-number)<0.0001){
      rounded+=1;
      new_coord=rounded*gridStep;
    }
  }else if(min_or_max==0){ // min value
    if(new_coord>number || fabs(new_coord-number)<0.0001){
      rounded-=1;
      new_coord=rounded*gridStep;
    }
  }
  
  return(new_coord);
}

//check if a coordinate is part of a grid size
int Grid::inGridRes(vertex& vrtx, float res){

  float x_diff,x_ref,y_diff,y_ref,z_diff,z_ref;
  
  x_diff=fabs(vrtx.x-min_x)/res;
  x_ref=round(x_diff);

  y_diff=fabs(vrtx.y-min_y)/res;
  y_ref=round(y_diff);

  z_diff=fabs(vrtx.z-min_z)/res;
  z_ref=round(z_diff);

  if((fabs(x_diff-x_ref)<0.0001) && (fabs(y_diff-y_ref)<0.0001) && (fabs(z_diff-z_ref)<0.0001)){
    return(1);
  }else{
    return(0);
  }
}

float dist_3d(float x1, float y1, float z1, float x2, float y2, float z2){
  float dist;
  double xdiff,ydiff,zdiff;

  xdiff=fabs((x2-x1));
  ydiff=fabs((y2-y1));
  zdiff=fabs((z2-z1));
  dist=sqrt((xdiff*xdiff)+(ydiff*ydiff)+(zdiff*zdiff));
  return(dist);
}

void Grid::smooth(){
  map<int,vertex>::iterator git;
  cout<<"Smoothing..."<<endl;
  for(git=this->GRID.begin(); git!=this->GRID.end(); git++){
    map<int,vertex>::iterator it;
    int id;
    vertex& m=git->second;
    if(m.p==1){ continue; }
    if(m.grid[0]!=1){ continue; }
    int xi=(int)((m.x-min_x)/gridStep)+1;
    int yi=(int)((m.y-min_y)/gridStep)+1;
    int zi=(int)((m.z-min_z)/gridStep)+1;


    for(int tx=xi-smoothDist; tx<=xi+smoothDist; tx++){
      id=generateID(width,height,tx,yi,zi);
      it = this->GRID.find(id);
      if(it != this->GRID.end()){
        if(this->GRID[id].p==1) continue;
        for(int j=0; j<6; j++){
          if(this->GRID[id].ints[j]==1){ m.ints[j]=1; }
        }
      }
    }
    for(int ty=yi-smoothDist; ty<=yi+smoothDist; ty++){
      id=generateID(width,height,xi,ty,zi);
      it = this->GRID.find(id);
      if(it != this->GRID.end()){
        if(this->GRID[id].p==1) continue;
        for(int j=0; j<6; j++){
          if(this->GRID[id].ints[j]==1){ m.ints[j]=1; }
        }
      }
    }
    for(int tz=zi-smoothDist; tz<=zi+smoothDist; tz++){
      id=generateID(width,height,xi,yi,tz);
      it = this->GRID.find(id);
      if(it != this->GRID.end()){
        if(this->GRID[id].p==1) continue;
        for(int j=0; j<6; j++){
          if(this->GRID[id].ints[j]==1){ m.ints[j]=1; }
        }
      }
    }
    
    //Diagonals
    // cout<<endl<<"going diag 0"<<endl;
    int tx;
    int ty;
    int tz;

    tx=xi-smoothDist; ty=yi+smoothDist; tz=zi+smoothDist;
    for(int s=0; s<6; s++){
      tx++; ty--; tz--;
      id=generateID(width,height,tx,ty,tz);
      it = this->GRID.find(id);
      if(it != this->GRID.end()){
        if(this->GRID[id].p==1) continue;
        for(int j=0; j<6; j++){
          if(this->GRID[id].ints[j]==1){ m.ints[j]=1; }
        }
      }
    }

    tx=xi-smoothDist; ty=yi+smoothDist; tz=zi-smoothDist;
    for(int s=0; s<6; s++){
      tx++; ty--; tz++;
      id=generateID(width,height,tx,ty,tz);
      it = this->GRID.find(id);
      if(it != this->GRID.end()){
        if(this->GRID[id].p==1) continue;
        for(int j=0; j<6; j++){
          if(this->GRID[id].ints[j]==1){ m.ints[j]=1; }
        }
      }
    }

    tx=xi+smoothDist; ty=yi+smoothDist; tz=zi-smoothDist;
    for(int s=0; s<6; s++){
      tx--; ty--; tz++;
      id=generateID(width,height,tx,ty,tz);
      it = this->GRID.find(id);
      if(it != this->GRID.end()){
        if(this->GRID[id].p==1) continue;
        for(int j=0; j<6; j++){
          if(this->GRID[id].ints[j]==1){ m.ints[j]=1; }
        }
      }
    }

    tx=xi-smoothDist; ty=yi-smoothDist; tz=zi-smoothDist;
    for(int s=0; s<6; s++){
      tx++; ty++; tz++;
      id=generateID(width,height,tx,ty,tz);
      it = this->GRID.find(id);
      if(it != this->GRID.end()){
        if(this->GRID[id].p==1) continue;
        for(int j=0; j<6; j++){
          if(this->GRID[id].ints[j]==1){ m.ints[j]=1; }
        }
      }
    }
  }
}

int is_coord_in_cube(float x, float y, float z, float center_x, float center_y, float center_z, float edge){
  float x_max,x_min,y_max,y_min,z_max,z_min;
  x_max=center_x+edge;
  x_min=center_x-edge;
  y_max=center_y+edge;
  y_min=center_y-edge;
  z_max=center_z+edge;
  z_min=center_z-edge;
  
  if( ((x<x_max) || (fabs((double)(x-x_max))<0.0001)) && ((x>x_min) || (fabs((double)(x-x_min))<0.0001))){ //if x is in the cube check y
    if( ((y<y_max) || (fabs((double)(y-y_max))<0.0001)) && ((y>y_min) || (fabs((double)(y-y_min))<0.0001))){ //if yes check z
      if( ((z<z_max) || (fabs((double)(z-z_max))<0.0001)) && ((z>z_min) || (fabs((double)(z-z_min))<0.0001))){ //if yes, ITS IN!!!
        return(1);
      }else{return(0);}
    }else{return(0);}
  }else{return(0);}
}

void Grid::writeMif(vector<atom>& prot, vector<atom>& lig){ //Write the .mif output file

  map<int,vertex>::iterator it;
  FILE* fpNew;
  int probe;
  float d=0.0;
  int bsFlag=0;
  int i,j;

  //Print header of output file
  string na="NA";
  string mifFileNew=outBase + tag + ".mif";
  fpNew = fopen(mifFileNew.c_str(),"w");
  fprintf(fpNew,"#cmd %s\n",cmdLine.c_str());
  fprintf(fpNew,"#proteinFile %s\n",proteinFile.c_str());
  fprintf(fpNew,"#cleftFile %s\n",cleftFile.c_str());
  if(resnumc.compare("")!=0){
    fprintf(fpNew, "#resnumc %s\n",resnumc.c_str());
  }else{
    fprintf(fpNew, "#resnumc none\n");
  }
  fprintf(fpNew, "#gridLigMinDist %5.2f\n",sqrt(gridLigDist));
  fprintf(fpNew,"#chain %s\n",chain.c_str());
  fprintf(fpNew,"#ff %s\n",ff.c_str());
  fprintf(fpNew,"#nbOfProbes %2d\n",nbOfProbes);
  for(i=0; i<nbOfProbes; i++){
    fprintf(fpNew,"#probe[%d] %4s\n",i,probes[i].c_str());
  }
  fprintf(fpNew,"#zip %d\n",zip);
  fprintf(fpNew,"#stepsize %4.2f\n",stepsize);
  fprintf(fpNew,"#gridStep %4.2f\n",gridStep);
  fprintf(fpNew,"#angThreshold %6.2f\n",angThresh);
  fprintf(fpNew,"#atom_probe_distance_threshold %7.4f\n",atmPbMaxDist);
  fprintf(fpNew,"#protein_grid_distance %7.4f to %7.4f\n",sqrt(minGridDist),sqrt(maxGridDist));
  fprintf(fpNew,"#grid_width %d\n",width);
  fprintf(fpNew,"#grid_height %d\n",height);
  fprintf(fpNew,"#gs %d %d %d %d\n",vrtx200,vrtx150,vrtx100,vrtx050);
  fprintf(fpNew,"#ss %d %d %d %d\n",ss[0],ss[1],ss[2],ss[3]);
  fprintf(fpNew,"#ssm %d %d %d %d\n",ssm[0],ssm[1],ssm[2],ssm[3]);

  cout<<endl<< "Writing Mif File"<<endl;
  for(it=GRID.begin();it!=GRID.end();it++){
    if(it->second.p==1){
      if(it->second.grid[2]==1){
        fprintf(fpNew, "#PG %7.2f %7.2f %7.2f\n",it->second.x,it->second.y,it->second.z);
      }
    }else{
      if((zip!=-1 && it->second.grid[zip]==1) || zip==-1){
        fprintf(fpNew, "%7.2f %7.2f %7.2f",it->second.x,it->second.y,it->second.z);
        for(probe=0; probe<nbOfProbes; probe++){
          fprintf(fpNew, " %1d",it->second.ints[probe]);
          if(it->second.nrgs[probe]==0 || abs(it->second.nrgs[probe]-0.0)<0.00001){
            fprintf(fpNew, " 0");
          }else{ fprintf(fpNew, " %6.3f",it->second.nrgs[probe]); }
          if(it->second.angles[probe]==0 || abs(it->second.angles[probe]-0.0)<0.00001){
            fprintf(fpNew, " 0");
          }else{ fprintf(fpNew, " %6.3f",it->second.angles[probe]); }
        }
        for(i=0; i<4; i++){
          fprintf(fpNew, " %1d",it->second.grid[i]);
        }
        // for(i=0; i<aa.size(); i++){
        //   if(it->second.env[aa[i]]>5.0){
        //     fprintf(fpNew, " %d",0);
        //   }else{
        //     fprintf(fpNew, " %d",(int)round(it->second.env[aa[i]]));
        //   }
        // }
        fprintf(fpNew, " %d\n",it->second.bu);
      }
    }
  }

  cout<<endl<< "Writing Pseudo"<<endl;
  for(int i=0; i<pseudoList.size(); i++){
    fprintf(fpNew,"#PSEUDO %s %8.3f %8.3f %8.3f\n",pseudoList[i].type.c_str(),pseudoList[i].x,pseudoList[i].y,pseudoList[i].z);
  }

  if(printAtoms==1){
    cout<<endl<< "Writing atoms"<<endl;
    if(zip==1){
      for(i=0; i<lig.size(); i++){
        fprintf(fpNew,"#ATOM %3s %4d %4s %5d %s %s %8.3f %8.3f %8.3f %d %d\n",lig[i].resn.c_str(),lig[i].resnb,lig[i].atomn.c_str(),lig[i].atomnb,lig[i].chain.c_str(),lig[i].alt.c_str(),lig[i].x,lig[i].y,lig[i].z,lig[i].mif,lig[i].bs);
      }
    }else{
      int step=(int)(caT/gridStep);
      for(i=0; i<prot.size(); i++){
        if(prot[i].mif!=1) continue;
        int xi=(int)((prot[i].x-min_x)/gridStep)+1;
        int yi=(int)((prot[i].y-min_y)/gridStep)+1;
        int zi=(int)((prot[i].z-min_z)/gridStep)+1;

        int flag=0;
        // cout<<prot[i].resn<<" "<<prot[i].resnb<<" "<<prot[i].chain<<" "<<prot[i].atomn<<" "<<prot[i].bs<<endl;
        for(int tx=xi-step; tx<=xi+step; tx++){
          for(int ty=yi-step; ty<=yi+step; ty++){
            for(int tz=zi-step; tz<=zi+step; tz++){
              int tid=generateID(width,height,tx,ty,tz);
              it = GRID.find(tid);
              if(it != GRID.end()){
                if(it->second.p==1) continue;
                if(it->second.grid[0]!=1) continue;
                // cout<<"exist "<<it->second.x<<" "<<it->second.y<<" "<<it->second.z<<endl;
                float d=sqrt(dist_3d(prot[i].x,prot[i].y,prot[i].z,it->second.x,it->second.y,it->second.z));
                // cout<<d<<endl;
                // cout<<it->second.x<<" "<<it->second.y<<" "<<it->second.z<<" "<<it->second.grid[0]<<" "<<it->second.grid[1]<<" "<<it->second.grid[2]<<" "<<it->second.grid[3]<< endl;
                if(d<caT || fabs(d-caT)<0.001){
                  prot[i].bs=1;
                  // cout<<prot[i].resn<<" "<<prot[i].resnb<<" "<<prot[i].chain<<endl;
                  flag=1;
                  break;
                }
              }
              if(flag==1) break;
            }
            if(flag==1) break;
          }
          if(flag==1) break;
        }
        // cout<<prot[i].resn<<" "<<prot[i].resnb<<" "<<prot[i].chain<<" "<<prot[i].atomn<<" "<<prot[i].bs<<endl;
      }  
      for(j=0; j<prot.size(); j++){
        fprintf(fpNew,"#ATOM %3s %4d %4s %5d %s %8.3f %8.3f %8.3f %d %d\n",prot[j].resn.c_str(),prot[j].resnb,prot[j].atomn.c_str(),prot[j].atomnb,prot[j].chain.c_str(),prot[j].x,prot[j].y,prot[j].z,prot[j].mif,prot[j].bs);
      }
    }   
  }
  fclose(fpNew);
}

void stripSpace(string &str) {
  for (int i=0;i<str.length();i++) 
    if (str[i]==' ') {
    str.erase(i,1);
    i--;
  }
}

void getAtomRef(){
  string fn=basePath + "/forcefield_files/"+ff+"/atoms";
  ifstream infile(fn.c_str());
  string s;
  string res;
  string atom;
  string type;

  while(getline(infile,s)){
    istringstream iss(s);
    vector<string> tokens;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));
    atomTypes[tokens[0] + "_" + tokens[1]]=tokens[2];

    angRef naf;
    naf.r1="";
    naf.r2="";
    if(tokens.size()==5){
      naf.r1=tokens[3];
      naf.rDir=atoi(tokens[4].c_str());
    }else if(tokens.size()==7){
      naf.r1=tokens[3];
      naf.r2=tokens[4];
      naf.rDir=atoi(tokens[5].c_str());
      naf.ring=atoi(tokens[6].c_str());
    }
    atomRef.insert(pair<string,angRef>(tokens[0] + "_" + tokens[1],naf));
  }
}

void getPseudoC(){
  string fn=basePath + "/forcefield_files/"+ff+"/pseudocenters";
  ifstream infile(fn.c_str());
  string line;
  string res;
  string atom;
  string pseudo;

  while(getline(infile,line)){
    stringstream test(line);
    test >> res >> atom >> pseudo;
    pseudoC[res+"_"+atom]=pseudo;
  }
}

void getEpsilons(){
  string fn;
  if(matrixF.compare("")!=0){
    fn=matrixF;
  }else{
    fn=basePath + "/forcefield_files/"+ff+"/epsilons";
  }
  ifstream infile(fn.c_str());
  string s;
  int row=0;
  vector<string> coln;

  while(getline(infile,s)){
    istringstream iss(s);
    vector<string> tokens;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));
    if(row==0){
      for(int i=0; i<tokens.size(); i++){
        coln.push_back(tokens[i]);
      }
    }else{
      string rown=tokens[0];
      for(int i=1; i<tokens.size(); i++){
        eps[coln[i-1]+"_"+rown]=atoi(tokens[i].c_str());
        eps[rown+"_"+coln[i-1]]=atoi(tokens[i].c_str());
      }
    }
    row++;
  }
}

void readLigFile(){
  ifstream infile(ligFile.c_str());
  string s;
  while(getline(infile,s)){
    istringstream iss(s);
    vector<string> tokens;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));
    // ATP 544 PG P.3
    // ATP 544 O1G O.2
    // ATP 544 O2G O.3
    ligAt[tokens[2]]=tokens[3];
    // cout<<tokens[0]<<" "<<tokens[1]<<" "<<tokens[2]<<" "<<ligAt[tokens[0]+"_"+tokens[1]+"_"+tokens[2]]<<endl;
  }
}

void getAtomTypes(){
  string fn=basePath + "/forcefield_files/"+ff+"/atomTypes";
  ifstream infile(fn.c_str());
  string s;
  int row=0;

  while(getline(infile,s)){
    // cout<<s<<endl;
    istringstream iss(s);
    vector<string> tokens;
    copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));
    row++;
    if(row==1) continue;
    don[tokens[1]]=atoi(tokens[2].c_str());
    acc[tokens[1]]=atoi(tokens[3].c_str());
    arm[tokens[1]]=atoi(tokens[4].c_str());
    chr[tokens[1]]=atoi(tokens[5].c_str());
    hyd[tokens[1]]=atoi(tokens[6].c_str());
    // cout<<tokens[1]<<" "<<don[tokens[1]]<<" "<<acc[tokens[1]]<<" "<<arm[tokens[1]]<<" "<<chr[tokens[1]]<<" "<<hyd[tokens[1]]<<endl;

  }
}

void getProbes(){
  string fn;
  if(probesF.compare("")!=0){
    fn=probesF;
  }else{
    fn=basePath + "/forcefield_files/"+ff+"/probes";
  }
  ifstream infile(fn.c_str());
  string line;
  string pbn;
  float min;
  float max;
  float thresh;

  while(getline(infile,line)){
    stringstream test(line);
    test >> pbn >> min >> max >> thresh;
    minD[pbn]=min;
    maxD[pbn]=max;
    nrgT[pbn]=thresh;
    // cout<< pbn <<" "<<minD[pbn]<<" "<<maxD[pbn]<<" "<<nrgT[pbn]<<endl;
    probes.push_back(pbn);
    idAt[nbOfProbes]=pbn;
    nbOfProbes++;
  }
}

void getaa(){
  string fn=basePath + "/forcefield_files/"+ff+"/aa";
  ifstream infile(fn.c_str());
  string line;
  string taa;

  while(getline(infile,line)){
    stringstream test(line);
    test >> taa;
    aa.push_back(taa);
  }
}

void getEnv(map<int,vertex>& grid, vector<atom>& prot, vector<int>& vrtxList){
  cout<<endl<< "Getting Environment"<< endl;
  for(int i=0; i<vrtxList.size(); i++){
    map<int,vertex>::iterator it=grid.find(vrtxList[i]);
    vertex& m=it->second;
    
    if(m.p==1){ continue; }
    if(m.grid[0]!=1 && m.grid[1]!=1 && m.grid[2]!=1 && m.grid[3]==1) continue;

    for(int j=0; j<prot.size(); j++){ //Iterate each atom for this probe at this grid intersection
      float dist=dist_3d(prot[j].x,prot[j].y,prot[j].z,m.x,m.y,m.z);
      if(dist<m.env[prot[j].resn]) m.env[prot[j].resn]=dist;
    }
  }
}

void getMif(map<int,vertex>& grid, vector<atom>& prot, vector<int>& vrtxList){
  int j;
  cout<<endl<< "Searching for potential interactions at each grid vertex"<< endl;

  for(int i=0; i<vrtxList.size(); i++){
    map<int,vertex>::iterator it=grid.find(vrtxList[i]);
    vertex& m=it->second;

    if(m.p==1){ continue; } //Skip if this vertex is a protein vertex

    if(m.grid[0]==1 || m.grid[1]==1 || m.grid[2]==1 || m.grid[3]==1){ //If its in one of the grid resolution (not sure if this is useful anymore)
    
      if(printDetails==1){ cout<<endl<<"Vertex id: "<< m.id <<" "<<m.x<<" "<<m.y<<" "<<m.z<<endl; }

      int flag=0;
      for(int probe=0; probe<nbOfProbes; probe++){ //Iterate each probe

        float enrg_sum=0.00;
        int countAtms=0;
        float closest=10000;

        if(printDetails==1){ cout<<endl<<"### PROBE "<<probe<<" ###"<<endl; }

        //Iterate each atom for this probe and sum the potential energy
        for(int j=0; j<prot.size(); j++){
          if(prot.at(j).mif!=1) continue;
          //Calculate and sum the potential energy
          enrg_sum+=calcNrg(m,prot.at(j),probe,countAtms,closest);
        }

        if(countAtms>0){
          m.nrgs[probe]=enrg_sum;
          if(enrg_sum<nrgT[probes[probe]] || (fabs(enrg_sum-nrgT[probes[probe]]))<0.001){
            m.ints[probe]=1;
            for(int gi=0; gi<4; gi++){
              if(m.grid[gi]==1){
                ssm[gi]++;
              }
            }
            flag=1;            
          }
          if(printDetails==1){ cout<<m.grid[0]<<" "<<m.grid[1]<<" "<<m.grid[2]<<" "<<m.grid[3]<<" probe "<<probe<<" nrg "<< enrg_sum<<" int "<<m.ints[probe]<<endl; }
        }else{
          if(printDetails==1){ cout<<"NRG "<<probe<<" no atoms"<<endl; }
        }
      }

      //Increment search space
      if(flag==1){
        for(int gi=0; gi<4; gi++){
          if(m.grid[gi]==1){
            ss[gi]++;
          }
        }
      }
    }
  }
}

double calcNrg(vertex& vrtx, atom& atm, int pbId, int& count_atoms, float& closest){
  float dist,alpha,epsilon,angle;
  int atomAtId,pbAtId;
  double energy;
  energy=0.0;
  alpha=1.0;
  angle=1.0;
  float rDist,rpDist;

  //Get distance between probe (vertex) and atom
  dist=dist_3d(atm.x,atm.y,atm.z,vrtx.x,vrtx.y,vrtx.z);

  //If atom is too far or too close to the probe specific limits or the general atmPbMaxDist limit return 0
  if(dist > maxD[probes[pbId]] || dist < minD[probes[pbId]] || dist > atmPbMaxDist) return(energy);

  string at=atomTypes[atm.resn+"_"+atm.atomn];
  string pat=probes[pbId];
  epsilon=eps[at+"_"+pat];

  if(printDetails==1){ cout<<atm.resn<<" "<<atm.resnb<<" "<<atm.atomn<<" "<<at<<" dist "<<dist<<endl; }

  //If its an acceptor or donor probe-atom interaction (defined in the input file)
  if(((acc[at]==1 && don[pat]==1) || (acc[pat]==1 && don[at]==1)) && atm.dir==1){

    //Calculate distances and angle using reference atoms
    rDist=dist_3d(atm.x,atm.y,atm.z,atm.xr,atm.yr,atm.zr);
    rpDist=dist_3d(vrtx.x,vrtx.y,vrtx.z,atm.xr,atm.yr,atm.zr);
    angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
    angle=acos(angle)* 180.0 / PI;

    if(atm.rDir==0){ angle=180.00-angle; }
    if(printDetails==1){ cout<<"Hbond Angle "<< angle <<endl; }

    //Check if angle is above threshold
    if(angle>angThresh){//The angle is too big
      if(printDetails==1){ cout<<"Angle over threshold"<<endl; }
      return(energy);
    }else{ //The angle is good
      if(printDetails==1){ cout<<"Hbond Angle "<< angle <<" threshold "<<angThresh<<endl; }
    }

    //If this is the closest donnor or acceptor atom to this probe store the hydrogen bond angle (that will be outputed in .mif file)
    if(dist<closest){
      closest=dist;
      vrtx.angles[pbId]=angle;
      if(printDetails==1){ cout<<" closest "<<dist<<" new angle "<<vrtx.angles[pbId]<<endl; }
    }

  }else if(arm[at]==1 && arm[pat]==1 && atm.dir==1){ //If its an aromatic interaction

    //Calculate distances and angle using reference atoms
    if(atm.dir==1){
      rDist=dist_3d(atm.x,atm.y,atm.z,atm.xr,atm.yr,atm.zr);
      rpDist=dist_3d(vrtx.x,vrtx.y,vrtx.z,atm.xr,atm.yr,atm.zr);
      angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
      angle=acos(angle)* 180.0 / PI;
      if(printDetails==1){ cout<<"Aromatic Angle "<< angle<<endl; }
      if(angle>90 || fabs(180.00-angle) < 0.001) angle=180-angle; //IF angle is above the threshold (which is 90 degrees by default, the equivalent of no angle threshold)
      // if(angle>25.00 && angle<65.00){
      //   if(printDetails==1){ cout <<"angle between 35 and 55"<<endl; }
        // return(energy);
      // }

      //If this is the closest aromatic atom to this probe store the hydrogen bond angle (that will be outputed in .mif file)
      if(dist<closest){
        closest=dist;
        vrtx.angles[pbId]=angle;
        if(printDetails==1){ cout<<" closest "<<dist<<" new angle "<<vrtx.angles[pbId]<<endl; }
      }
    }else{
      if(printDetails==1){ cout<<"Both aromatic but no angle to calculate."<<endl; }
    }
  }else if(chr[at]==1 && chr[pat]==1){ //charged interaction
    if(printDetails==1){ cout<<"Charged couple "<<endl; }
  }else if(hyd[pat]==1 && hyd[at]==1){//if its a hydrophoic probe
    if(printDetails==1){ cout << "Hydrophobic "<<endl; }
  }
  count_atoms++;
  energy=(epsilon)*(exp(-1.0*dist*alpha));
  if(printDetails==1){ cout<<"Epsilon "<<epsilon<<" NRG "<<energy<<endl; }
  return(energy);
}

void getPseudo(map<int,vertex>& grid, vector<atom>& prot, vector<int>& vrtxList){
  cout<<endl<< "Projecting pseudocenters..."<< endl;

  for(int j=0; j<prot.size(); j++){
    if(prot[j].mif!=1) continue;

    // cout<< prot[j].resn + " " + prot[j].atomn + " " << prot[j].resnb<< " mif "<< prot[j].mif;

    map<string,string>::iterator it = pseudoC.find(prot[j].resn+"_"+prot[j].atomn);
    if(it != pseudoC.end()){

      string pseudo=it->second;
      int angleType=0;
      if(pseudo.compare("doa")==0 || pseudo.compare("don")==0 || pseudo.compare("acc")==0){
        angleType=1;
      }else if(pseudo.compare("arm")==0){
        // angleType=2;
      }
      // cout<<" "+pseudo<<" angletype "<<angleType<<endl;
    
      int bestID=-1;
      float bestD=10.0;
      for(int i=0; i<vrtxList.size(); i++){
        vertex& m=grid.at(vrtxList.at(i));
        if(m.bu<bul) continue;
        if(m.grid[2]!=1) continue;

        float angle=0.0;
        float dist=dist_3d(prot[j].x,prot[j].y,prot[j].z,m.x,m.y,m.z);
        
        if(dist>4.0) continue;
        if(angleType==1){
          float rDist=dist_3d(prot[j].x,prot[j].y,prot[j].z,prot[j].xr,prot[j].yr,prot[j].zr);
          float rpDist=dist_3d(m.x,m.y,m.z,prot[j].xr,prot[j].yr,prot[j].zr);
          angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
          angle=acos(angle)* 180.0 / PI;
          if(prot[j].rDir==0){ angle=180.00-angle; }
          if(angle > 60.00) continue;
          // cout<< " angle "<<angle<<" dist "<<dist<<endl;
        }else if(angleType==2){ //aromatic angle
          // float rDist=dist_3d(prot[j].x,prot[j].y,prot[j].z,prot[j].xr,prot[j].yr,prot[j].zr);
          // float rpDist=dist_3d(m.x,m.y,m.z,prot[j].xr,prot[j].yr,prot[j].zr);
          // angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
          // angle=acos(angle)* 180.0 / PI;
          // if(angle>90 || fabs(180.00-angle) < 0.001){
          //   angle=180-angle;
          // }
          // cout<< " angle "<<angle<<" dist "<<dist<<endl;
        }
        if(dist<bestD){
          bestD=dist;
          bestID=i;
        }
      }
      if(bestID!=-1){
        pseudovrtx npv;
        npv.dist=bestD;
        npv.type=pseudo;
        npv.x=grid.at(vrtxList.at(bestID)).x;
        npv.y=grid.at(vrtxList.at(bestID)).y;
        npv.z=grid.at(vrtxList.at(bestID)).z;
        pseudoList.push_back(npv);
      }
    }
  }
}

void getStats(map<int,vertex>& grid, vector<atoms>& lig,vector<int>& vrtxList){
  string newFile=statsF + tag;
  FILE* fpNew;
  // fpNew = fopen(newFile.c_str(),"w");
  for(int probe=0; probe<nbOfProbes; probe++){ //Iterate each probe
    cout<<probe<<endl;
    int tp=0;
    int tn=0;
    vector<pbVrtx> pbV;
    for(int l=0; l<lig.size(); l++){
      if(lig[l].at.compare(idAt[probe])==0){
        for(int i=0; i<vrtxList.size(); i++){
          map<int,vertex>::iterator it=grid.find(vrtxList[i]);
          vertex& m=it->second;

          if(m.p==1) continue;
          if((zip!=-1 && it->second.grid[zip]==1) || zip==-1){
            float dist=dist_3d(lig[l].x,lig[l].y,lig[l].z,m.x,m.y,m.z);
            // fprintf(fpNew,"%8s %8.3f %11.4f\n",idAt[probe].c_str(),dist,m.nrgs[probe]);
            pbVrtx npv;
            npv.dist=dist;
            npv.nrg=m.nrgs[probe];
            if(dist<distpbV){ tp++; }else{ tn++; }
            pbV.push_back(npv);
            // cout<<"Vertex id: "<< m.id <<" "<<m.x<<" "<<m.y<<" "<<m.z<<endl;
            // cout<<"probe "<<probe<<" "<<idAt[probe]<<" "<<m.nrgs[probe]<<" dist "<<dist<<" "<<lig[l].at<<endl;
          }
        }
      }
    }
    if(pbV.size()!=0){
      sort(pbV.begin(), pbV.end(), compByNrg);
      int ctp=0;
      int ctn=0;
      int rocV[tn];
      rocV[0]=0;
      for(int i=0; i<pbV.size(); i++){
        cout<<idAt[probe]<<" "<<pbV.at(i).dist<<" "<<pbV.at(i).nrg;
        if(pbV.at(i).dist<distpbV){
          cout<< " TP";
          ctp++;
        }else{
          ctn++;
        }
        cout<<endl;
        rocV[ctn]=ctp;
      }
      float xs=1.0/ctn;
      float ys=1.0/ctp;
      float roc=0.0;
      // cout<<"ctp "<<ctp<<" ctn "<<ctn<<" xs "<<xs<<" ys "<<ys<<endl;
      for(int j=1; j<=tn; j++){
        roc+=(xs*((float)rocV[j-1])*ys)+(xs*(((float)rocV[j]*ys)-((float)rocV[j-1]*ys)))/2;
        // cout<<j<<" "<<rocV[j]<<endl;
      }
      cout<<"roc "<<roc<<" tp "<<ctp<<" ctn "<<ctn<<endl;
    }
  }
  // fclose(fpNew);
}

bool compByNrg(pbVrtx a, pbVrtx b)
{
    return a.nrg < b.nrg;
}

void getPairwise(){
  cout<<"Getting pariwise"<<endl;
  string empt ("");
  string line;
  ifstream infile(pairwiseF.c_str());
  while(getline(infile,line)){
    if(line.compare("")!=0){
      istringstream ss(line);
      istream_iterator<string> begin(ss), end;
      vector<string> vec(begin, end);
      string pdbF="";
      string cleftF="";
      string rnc="";
      string ligF="";
      for(int i=0; i<vec.size(); i+=2){
        if(vec[i].compare("-p")==0) pdbF=vec[i+1];
        if(vec[i].compare("-g")==0) cleftF=vec[i+1];
        if(vec[i].compare("-l")==0) rnc=vec[i+1];
        if(vec[i].compare("-lf")==0) ligF=vec[i+1];
      }
      pwRun npw;
      npw.pdbF=pdbF;
      npw.cleftF=cleftF;
      npw.rnc=rnc;
      npw.ligF=ligF;
      pw.push_back(npw);
    }
  }
}
