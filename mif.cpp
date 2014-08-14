//By Matthieu Chartier

#include "mif.h"
#include "./forcefield_files/getAtomId_FAINT.h"

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

  if(readCmdLine(argc, argv)==24){ return(0); }

  NEW(atPair,atPairStruc);
  NEW(atomSingle,atomSingleStruc);
  NEW(atSingle,atSingleStruc);

  storeAtomSingle_FAINT(atomSingle);
  storeAtSingle_FAINT(atSingle);
  storeAtPair_FAINT(atPair);
  getProbes_FAINT();

  //Get prefix
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

  Protein protein=Protein(proteinFile);

	Grid grid=Grid(cleftFile,protein);
  cout<< "Grid has " <<grid.GRID.size()<< " vertexes."<<endl;
  cout<< "Vertex list size: "<< grid.vrtxIdList.size()<<endl;

  get_enrg(grid.GRID,protein.PROTEINOBJ,grid.vrtxIdList);
  //if(stepsize < 2.0){ grid.smooth(); }
  grid.writeMif(protein.PROTEIN);

	return(0);
}


//-------------------------------END OF MAIN-----------------------------
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/

Protein::~Protein(void){

}

Protein::Protein(string filename){
  readPDB(filename);
  createProteinObject();
}

int readCmdLine(int argc, char **argv){
  stringstream usage;
  int nb_arg;

  for(nb_arg=0; nb_arg<argc; nb_arg++){
    cmdLine = cmdLine + argv[nb_arg];
    if(nb_arg<argc-1){
      cmdLine = cmdLine + " ";
    }
  }
    
  usage << "\n!---   IsoMIF   ---!\nWelcome.Bienvenue\n";
  usage << "\nObligatory Arguments\n";
  usage << "-p <filename>        : \t Protein (PDB format)\n";
  usage << "-g <filename>        : \t Cleft of protein (GetCleft sphere format)\n";
  usage << "\nOptional Arguments\n";
  usage << "-s                   : \t grid spacing 0.25, 0.5, 1.0 or 2.0 (default: 0.25Angstrom)\n";
  usage << "-d <distance in Ang.>: \t atom-probe distance threshold. (default: "<< atmPbMaxDist <<" Angstrom)\n";
  usage << "-c <chain>           : \t which Chain to consider for energy calculations\n";
  usage << "-w <distance in Ang.>: \t minimum distance of the probes to the atoms (default: "<< minGridDist <<" Angstrom)\n";
  usage << "-v <distance in Ang.>: \t maximum distance of the probes to the atoms (default: "<< maxGridDist <<" Angstrom)\n";
  usage << "-pr                  : \t print details of energy calculations for each probe-atom pair\n";
  usage << "-o                   : \t output directory\n";
  usage << "-t                   : \t filename (tag) [pdb name by default]\n";
  usage << "-l                   : \t RESNUMC of the ligand from which to crop the grid\n";
  usage << "-r                   : \t maximum distance between the grid and the ligand\n";
  usage << "-x                   : \t do not write atoms in mif file\n";
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
    if(strcmp(argv[nb_arg],"-c")==0){
      chain=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-l")==0){
      resnumc=string(argv[nb_arg+1]);
    }
    if(strcmp(argv[nb_arg],"-v")==0){
      sscanf(argv[nb_arg+1], "%f", &maxGridDist);
    }
    if(strcmp(argv[nb_arg],"-w")==0){
      sscanf(argv[nb_arg+1], "%f", &minGridDist);
    }
    if(strcmp(argv[nb_arg],"-r")==0){
      sscanf(argv[nb_arg+1], "%f", &gridLigDist);
    }
    if(strcmp(argv[nb_arg],"-d")==0){
      sscanf(argv[nb_arg+1], "%f", &atmPbMaxDist);
    }
    if(strcmp(argv[nb_arg],"-s")==0){
      sscanf(argv[nb_arg+1], "%f", &stepsize);
    }
    if(strcmp(argv[nb_arg],"-pr")==0){
      printDetails=1;
    }
    if(strcmp(argv[nb_arg],"-t")==0){
      tag=string(argv[nb_arg+1]);
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

  cout<<endl;
  cout<< "ProteinFile: "<< proteinFile <<endl;
  cout<< "CleftFile: "<< cleftFile <<endl;
  cout<< "Outbase: "<< outBase <<endl;
  cout<< "Chain: "<< chain <<endl;
  cout<< "Tag: " << tag << endl;
  cout<< "Stepsize: "<<stepsize<<endl;
  cout<< "MinGridDist: "<< minGridDist <<endl;
  cout<< "MaxGridDist: "<< maxGridDist <<endl;
  cout<< "AtmPbMaxDist: "<< atmPbMaxDist <<endl;
  cout<< "RESNUMC: "<< resnumc <<endl;
  cout<< "gridLigDist: "<< gridLigDist <<endl;
  gridLigDist=gridLigDist*gridLigDist;
  maxGridDist=maxGridDist*maxGridDist;
  minGridDist=minGridDist*minGridDist;

  return(0);
}

void Protein::readPDB(string filename){
  ifstream ifs;
  ofstream ofs;
  string line;
  string fields[13];
  string thisresnumc;
  float x,y,z;
  int atomnb;
  int resnb;
  atom* atm=NULL;
  string outFileName=outBase + tag + "_cpy.pdb";
  size_t found;
  ifs.open(filename.c_str());
  ofs.open(outFileName.c_str());

  if(!ifs.is_open()){ 
    cout << "could not read "<< filename << endl;
  }
  cout<<"Reading PDB"<<endl;
  while(ifs.good()){
    getline(ifs,line);

    if(line.compare(0,3,"END") == 0){ break; }

    //Store ligand if necessary using resnumc
    if(line.compare(0,6,"HETATM") == 0 && resnumc.compare("")!=0){
      thisresnumc = line.substr(17,3) + line.substr(22,4) + line.substr(21,1);
      stripSpace(thisresnumc);
      // cout<<endl<<thisresnumc;
      if(resnumc.compare(thisresnumc)==0){
        LIGAND.push_back(atof((line.substr(30,8).c_str())));
        LIGAND.push_back(atof((line.substr(38,8).c_str())));
        LIGAND.push_back(atof((line.substr(46,8).c_str())));
        // cout<< resnumc<< " to "<< thisresnumc << " "<< atof((line.substr(30,8).c_str()))<<" "<< atof((line.substr(38,8).c_str()))<<" "<<atof((line.substr(46,8).c_str()))<<endl;
      }
    }

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
      stripSpace(fields[6]);
      atomnb=atoi(fields[1].c_str());
      resnb=atoi(fields[6].c_str());
      x=atof(fields[8].c_str());
      y=atof(fields[9].c_str());
      z=atof(fields[10].c_str());
      // cout<<line<<endl;

      found = fields[2].find("H");

      if(resnumc.compare("")!=0){
        thisresnumc = line.substr(17,3) + line.substr(22,4) + line.substr(21,1);
        stripSpace(thisresnumc);
        if(resnumc.compare(thisresnumc)==0 && found==string::npos){
          LIGAND.push_back(atof((line.substr(30,8).c_str())));
          LIGAND.push_back(atof((line.substr(38,8).c_str())));
          LIGAND.push_back(atof((line.substr(46,8).c_str())));
          cout<< resnumc<< " to "<< thisresnumc << " "<< atof((line.substr(30,8).c_str()))<<" "<< atof((line.substr(38,8).c_str()))<<" "<<atof((line.substr(46,8).c_str()))<<endl;
        }
      }

      if((fields[3].compare("A")!=0) && (fields[3].compare("")!=0)){ continue; }
      if(fields[4].compare("HOH")==0){ continue; }

      atom atm;
      atm.x=x;
      atm.y=y;
      atm.z=z;
      atm.chain=fields[5];
      atm.atomnb=atomnb;
      atm.resnb=resnb;
      atm.atomn=fields[2];
      atm.resn=fields[4];
      atm.bs=0;
      atm.mif=1;

      if(found!=string::npos){
        atm.h=1;
      }else{
        atm.h=0;
      }

      if(getAtomId_FAINT(fields[4],fields[2])!=-1){
        atm.atomId=getAtomId_FAINT(fields[4],fields[2]);
      }else{
        atm.atomId=99999;
      }

      if(fields[5].compare(chain)!=0 && chain.compare("none")!=0){ atm.mif=0; }
      if(line.compare(0,6,"HETATM") == 0){ atm.mif=0; }

      PROTEIN.push_back(atm);
    }
  }
  ifs.close();
  ofs.close();
}

void Protein::createProteinObject(){
  int i,found,ring,needRef,rDir;
  float xr,yr,zr;
  string tatomn;
  string tatomn2;

  for(i=0; i<PROTEIN.size(); i++){
    if(PROTEIN[i].atomId!=99999 && PROTEIN[i].mif==1){
      ring=0;
      found=0;
      needRef=1;
      if(PROTEIN[i].resn.compare("ALA")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("ALA")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ARG")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("ARG")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ARG")==0 && PROTEIN[i].atomn.compare("NE")==0){
        tatomn="HE";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("ARG")==0 && PROTEIN[i].atomn.compare("NH1")==0){
        tatomn="CZ";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ARG")==0 && PROTEIN[i].atomn.compare("NH2")==0){
        tatomn="CZ";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ASN")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("ASN")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ASN")==0 && PROTEIN[i].atomn.compare("OD1")==0){
        tatomn="CG";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ASN")==0 && PROTEIN[i].atomn.compare("ND2")==0){
        tatomn="CG";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ASP")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("ASP")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ASP")==0 && PROTEIN[i].atomn.compare("OD1")==0){
        tatomn="CG";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ASP")==0 && PROTEIN[i].atomn.compare("OD2")==0){
        tatomn="CG";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("CYS")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("CYS")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("GLN")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("GLN")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("GLN")==0 && PROTEIN[i].atomn.compare("OE1")==0){
        tatomn="CD";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("GLN")==0 && PROTEIN[i].atomn.compare("NE2")==0){
        tatomn="CD";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("GLU")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("GLU")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("GLU")==0 && PROTEIN[i].atomn.compare("OE1")==0){
        tatomn="CD";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("GLU")==0 && PROTEIN[i].atomn.compare("OE2")==0){
        tatomn="CD";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("GLY")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("GLY")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("HIS")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("HIS")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("HIS")==0 && PROTEIN[i].atomn.compare("ND1")==0){
        tatomn="CG";
        tatomn2="CE1";
        ring=2;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("HIS")==0 && PROTEIN[i].atomn.compare("NE2")==0){
        tatomn="CE1";
        tatomn2="CD2";
        ring=2;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("HIS")==0 && (PROTEIN[i].atomn.compare("CD2")==0 || PROTEIN[i].atomn.compare("CE1")==0 || PROTEIN[i].atomn.compare("CG")==0)){
        tatomn="ND1";
        tatomn2="NE2";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("ILE")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("ILE")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("LEU")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("LEU")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("LYS")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("LYS")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("LYS")==0 && PROTEIN[i].atomn.compare("NZ")==0){
        tatomn="CE";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("MET")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("MET")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("PHE")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("PHE")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("PHE")==0 && PROTEIN[i].atomn.compare("CG")==0){
        tatomn="CD1";
        tatomn2="CD2";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("PHE")==0 && PROTEIN[i].atomn.compare("CD1")==0){
        tatomn="CG";
        tatomn2="CD2";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("PHE")==0 && (PROTEIN[i].atomn.compare("CD2")==0 || PROTEIN[i].atomn.compare("CE1")==0 || PROTEIN[i].atomn.compare("CE2")==0 || PROTEIN[i].atomn.compare("CZ")==0)){
        tatomn="CG";
        tatomn2="CD1";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("PRO")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("SER")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("SER")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("SER")==0 && PROTEIN[i].atomn.compare("OG")==0){
        tatomn="CB";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("THR")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("THR")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("THR")==0 && PROTEIN[i].atomn.compare("OG1")==0){
        tatomn="CB";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TRP")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("TRP")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TRP")==0 && PROTEIN[i].atomn.compare("CG")==0){
        tatomn="CD1";
        tatomn2="CD2";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TRP")==0 && PROTEIN[i].atomn.compare("CD1")==0){
        tatomn="CG";
        tatomn2="CD2";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TRP")==0 && (PROTEIN[i].atomn.compare("CD2")==0 || PROTEIN[i].atomn.compare("CE2")==0 || PROTEIN[i].atomn.compare("CE3")==0 || PROTEIN[i].atomn.compare("CZ2")==0 || PROTEIN[i].atomn.compare("CZ3")==0 || PROTEIN[i].atomn.compare("CH2")==0)){
        tatomn="CG";
        tatomn2="CD1";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TRP")==0 && PROTEIN[i].atomn.compare("NE1")==0){
        tatomn="CE2";
        tatomn2="CD1";
        ring=2;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TYR")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("TYR")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TYR")==0 && PROTEIN[i].atomn.compare("CG")==0){
        tatomn="CD1";
        tatomn2="CD2";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TYR")==0 && PROTEIN[i].atomn.compare("CD1")==0){
        tatomn="CG";
        tatomn2="CD2";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TYR")==0 && (PROTEIN[i].atomn.compare("CD2")==0 || PROTEIN[i].atomn.compare("CE1")==0 || PROTEIN[i].atomn.compare("CE2")==0 || PROTEIN[i].atomn.compare("CZ")==0)){
        tatomn="CG";
        tatomn2="CD1";
        ring=1;
        rDir=0;
      }else if(PROTEIN[i].resn.compare("TYR")==0 && PROTEIN[i].atomn.compare("OH")==0){
        tatomn="CZ";
        rDir=0;
      }else if(PROTEIN[i].resn.compare("VAL")==0 && PROTEIN[i].atomn.compare("N")==0){
        tatomn="H";
        rDir=1;
      }else if(PROTEIN[i].resn.compare("VAL")==0 && PROTEIN[i].atomn.compare("O")==0){
        tatomn="C";
        rDir=0;
      }else{
        needRef=0;
      }

      if(needRef==1){
        if(getRefAtom(xr, yr, zr, PROTEIN[i].resn, PROTEIN[i].resnb, tatomn, tatomn2, ring, PROTEIN[i].x, PROTEIN[i].y, PROTEIN[i].z, PROTEIN[i].chain)==1){
          atom atm;
          atm.x=PROTEIN[i].x;
          atm.y=PROTEIN[i].y;
          atm.z=PROTEIN[i].z;
          atm.atomId=PROTEIN[i].atomId;
          atm.chain=PROTEIN[i].chain;
          atm.atomnb=PROTEIN[i].atomnb;
          atm.resnb=PROTEIN[i].resnb;
          atm.atomn=PROTEIN[i].atomn;
          atm.resn=PROTEIN[i].resn;
          atm.xr=xr;
          atm.yr=yr;
          atm.zr=zr;
          atm.dir=1;
          atm.rDir=rDir;
          atm.bs=PROTEIN[i].bs;
          atm.mif=PROTEIN[i].mif;
          PROTEINOBJ.push_back(atm);
          // cout<< atm.resn << " " << atm.resnb << " " << atm.atomn << " " << atm.atomnb<<" "<<atm.atomId << " "<<atm.x << " "<<atm.y<< " "<<atm.z<< " |ref: "<<xr << " "<<yr<< " "<<zr<<" dir "<<atm.rDir<< endl;
        }
      }else{
        atom atm;
        atm.x=PROTEIN[i].x;
        atm.y=PROTEIN[i].y;
        atm.z=PROTEIN[i].z;
        atm.atomId=PROTEIN[i].atomId;
        atm.chain=PROTEIN[i].chain;
        atm.atomnb=PROTEIN[i].atomnb;
        atm.resnb=PROTEIN[i].resnb;
        atm.atomn=PROTEIN[i].atomn;
        atm.resn=PROTEIN[i].resn;
        atm.bs=PROTEIN[i].bs;
        atm.mif=PROTEIN[i].mif;
        atm.xr=0.0;
        atm.yr=0.0;
        atm.zr=0.0;
        atm.dir=0;
        PROTEINOBJ.push_back(atm);
      }
    }
  }
  cout<<endl<<"Protein has "<< PROTEINOBJ.size()<<" atoms"<<endl;
}

int Protein::getRefAtom(float& xr, float& yr, float& zr, string tresn, int tresnb, string tatomn, string tatomn2, int ring, float x, float y, float z, string chain){
  int i;
  int found=0;
  int foundr1=0;
  int foundr2=0;
  float xr1,yr1,zr1,xr2,yr2,zr2;
  float ax,ay,az,bx,by,bz;
  float dist,rDist,rpDist,angle;

  if(ring==0){
    for(i=0; i<PROTEIN.size(); i++){
      if(PROTEIN[i].resn.compare(tresn)==0 && PROTEIN[i].atomn.compare(tatomn)==0 && PROTEIN[i].resnb==tresnb && PROTEIN[i].chain.compare(chain)==0){
        found=1;
        xr=PROTEIN[i].x;
        yr=PROTEIN[i].y;
        zr=PROTEIN[i].z;
      }
      if(found==1){ break; }
    }
  }else{
    for(i=0; i<PROTEIN.size(); i++){
      if(PROTEIN[i].resn.compare(tresn)==0 && PROTEIN[i].atomn.compare(tatomn)==0 && PROTEIN[i].resnb==tresnb && PROTEIN[i].chain.compare(chain)==0){
        foundr1=1;
        xr1=PROTEIN[i].x;
        yr1=PROTEIN[i].y;
        zr1=PROTEIN[i].z;
      }
      if(PROTEIN[i].resn.compare(tresn)==0 && PROTEIN[i].atomn.compare(tatomn2)==0 && PROTEIN[i].resnb==tresnb && PROTEIN[i].chain.compare(chain)==0){
        foundr2=1;
        xr2=PROTEIN[i].x;
        yr2=PROTEIN[i].y;
        zr2=PROTEIN[i].z;
      }
      if(foundr1==1 && foundr2==1){
        break;
      }
    }
    if(foundr1==1 && foundr2==1){
      found=1;
      if(ring==1){
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
      }else if(ring==2){
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
  
  if(found==0){
    // if(ring==1){
    //   cout<<"Couldn't find ring vector of "<<tresn<<" "<<tresnb<<endl;
    // }else{
    //   cout<<"Couldn't find vector of "<<tresn<<" "<<tresnb<<" "<< tatomn<<endl;
    // }
  }else{
    if(ring==1){
      // cout<<endl<<"Looking for ring vector of "<<tresn<<" "<<tresnb<<endl;
      // cout<<"found "<<tatomn<<" "<<xr1<<" "<<yr1<<" "<<zr1<<endl;
      // cout<<"found "<<tatomn2<<" "<<xr2<<" "<<yr2<<" "<<zr2<<endl;
      // cout<<"base "<<x<<" "<<y<<" "<<z<<endl;
      // cout<<"ref "<< xr <<" "<< yr <<" "<< zr <<endl;
      // cout<<"angle "<<angle<<endl;
    }else if(ring==2){
      // cout<<endl<<"Looking for histidine vector of "<<tresn<<" "<<tresnb<<endl;
      // cout<<"found "<<tatomn<<" "<<xr1<<" "<<yr1<<" "<<zr1<<endl;
      // cout<<"found "<<tatomn2<<" "<<xr2<<" "<<yr2<<" "<<zr2<<endl;
      // cout<<"base "<<x<<" "<<y<<" "<<z<<endl;
      // cout<<"ref "<< xr <<" "<< yr <<" "<< zr <<endl;
      // cout<<"angle "<<angle<<endl;
    }else{
      // cout<<endl<<"Looking for atom vector of "<<tresn<<" "<<tresnb<<" "<< tatomn<<endl;
      // cout<<"found it at " << xr << " "<< yr << " " << zr << endl;
    }
  }

  return(found);
}

Grid::Grid(string filename, Protein& prot){
  min_x=1000.00;
  min_y=1000.00;
  min_z=1000.00;
  max_x=-1000.00;
  max_y=-1000.00;
  max_z=-1000.00;
  vrtx025=0;
  vrtx050=0;
  vrtx100=0;
  vrtx150=0;
  vrtx200=0;

  getMinMax(filename);
  readGetCleft(filename, prot.PROTEIN, prot.LIGAND);
}

Grid::~Grid(void){
}

int Grid::generateID(int w, int h, int x, int y, int z){
  int id;
  id=z * (w * h) + y * w + x;
  return(id);
}

void Grid::getMinMax(string filename){
  ifstream ifs;
  string line;
  string fields[4];
  float x,y,z,rad;
  float minx,miny,minz,maxx,maxy,maxz;

  ifs.open(filename.c_str());
  if(!ifs.is_open()){ 
    cout <<"could not read PDB file " << filename << endl;
  }

  while(ifs.good()){
    getline(ifs,line);
        
    if(line.compare(0,6,"ATOM  ") == 0){
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


      //Set GRID min/max X,Y and Z
      if(minx<min_x){ min_x=minx; }
      if(miny<min_y){ min_y=miny; }
      if(minz<min_z){ min_z=minz; }
      if(maxx>max_x){ max_x=maxx; }
      if(maxy>max_y){ max_y=maxy; }
      if(maxz>max_z){ max_z=maxz; }
    }
  }

  //Set GRID width and height
  width=(int)(((max_x-min_x)/stepsize)+1.0);
  height=(int)(((max_y-min_y)/stepsize)+1.0);

  cout <<endl<< "GRID min/max values: minx: "<< min_x<< " miny: "<< min_y<< " minz: "<<min_z<<" maxx: "<<max_x<<" maxy: "<<max_y<<" maxz: "<<max_z<<endl;
  cout<< "GRID Width "<<width<<" Height "<< height << endl;
}  

int Grid::readGetCleft(string filename, vector<atom>& protVec, vector<float>& ligVec){
  ifstream ifs;
  string line="";
  string fields[4];
  map<int,vertex>::iterator it;
  vector<atom>::iterator pit;
  vertex* vrtx=NULL;
  int i=0;
  int id;
  int uID=0;
  float x,y,z,rad,nx,ny,nz;
  float minx,miny,minz,maxx,maxy,maxz;
  float dist,minDist;

  ifs.open(filename.c_str());

  if(!ifs.is_open()){ 
    cout << "could not read "<< filename << endl;
  }

  while(ifs.good()){
    getline(ifs,line);
        
    if(line.compare(0,6,"ATOM  ") == 0 || line.compare(0,6,"HETATM") == 0){
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

      for(nx=minx; nx<=maxx; nx+=stepsize){
        for(ny=miny; ny<=maxy; ny+=stepsize){
          for(nz=minz; nz<=maxz; nz+=stepsize){

            //Generate grid vertex ID
            id=generateID(width,height,(int)((nx-min_x)/stepsize)+1,(int)((ny-min_y)/stepsize)+1,(int)((nz-min_z)/stepsize)+1);
            it = GRID.find(id);
            
            if(it != GRID.end()){
              continue;
            }

            //Skip if not within the sphere
            if(((abs(nx-x)*abs(nx-x))+(abs(ny-y)*abs(ny-y))+(abs(nz-z)*abs(nz-z)))>(rad*rad)){
              continue;
            }

            //Skip if too far from the ligand, if necessary
            if(resnumc.compare("")!=0){
              //cout<<"entering!"<<endl;
              minDist=10000.0;
              for(i=0; i<ligVec.size(); i+=3){
                dist=((abs(nx-ligVec.at(i))*abs(nx-ligVec.at(i)))+(abs(ny-ligVec.at(i+1))*abs(ny-ligVec.at(i+1)))+(abs(nz-ligVec.at(i+2))*abs(nz-ligVec.at(i+2))));
                //cout<< ligVec.at(i) << " "<< ligVec.at(i+1) << " "<< ligVec.at(i+2) << " grid: "<<nx<<" "<<ny<<" "<<nz<< " dist:"<< dist<< endl;
                if(dist<minDist){
                  minDist=dist;
                  //cout<< "New mindist: "<< minDist<<endl;
                }
              }

              if(minDist>gridLigDist){ //Skip to next grid intersection if too far from the ligand
                continue;
              }
            }

            //Skip if too far/close to the protein
            minDist=10000.0;
            for(pit=protVec.begin(); pit!=protVec.end(); ++pit){
              if((*pit).h==1) continue;
              if((*pit).mif==0) continue;
              dist=((abs(nx-(*pit).x)*abs(nx-(*pit).x))+(abs(ny-(*pit).y)*abs(ny-(*pit).y))+(abs(nz-(*pit).z)*abs(nz-(*pit).z)));
              if(dist<minDist){ minDist=dist; }
            }
            if(minDist<minGridDist || minDist > maxGridDist){
              continue;
            }
            
            if(it == GRID.end()){ //If this grid point doesnt exist
              vertex vrtx;
              vrtx.x=nx;
              vrtx.y=ny;
              vrtx.z=nz;

              vrtx.ints=new int[nbOfProbes];              
              if(vrtx.ints==NULL){
                printf("\n\nCan't malloc int**\nGoodbye.\n");
                return(24);
              }
              
              for(i=0; i<nbOfProbes; i++){ vrtx.ints[i]=0; }

              //cout<<"exist "<< vrtx.x << " "<<vrtx.y<<" "<<vrtx.z<<" "<<endl;

              //Print grid in appropriate file
              if(inGridRes(vrtx,2.0)==1){
                vrtx.grid[0]=1;
                vrtx200++;
              }else{ vrtx.grid[0]=0; }

              //Print grid in appropriate file
              if(inGridRes(vrtx,1.5)==1){
                vrtx.grid[1]=1;
                vrtx150++;
              }else{ vrtx.grid[1]=0; }
              
              if(inGridRes(vrtx,1.0)==1){
                vrtx.grid[2]=1;
                vrtx100++;
              }else{ vrtx.grid[2]=0; }
              
              if(inGridRes(vrtx,0.5)==1){
                vrtx.grid[3]=1;
                vrtx050++;               
              }else{ vrtx.grid[3]=0; }
              vrtx025++;

              vrtx.id=uID;
              uID++;
              GRID.insert(pair<int,vertex>(id,vrtx));
              vrtxIdList.push_back(id);
            }
          }
        }
      }
    }
  }
  ifs.close();

  cout<<"Grid volumes "<<vrtx200<<" "<<vrtx150<<" "<<vrtx100<<" "<<vrtx050<<endl;
  return(0);
}

float Grid::roundCoord(float number, int min_or_max){ 
  //min_or_max (1=max, 0=min)
  float quotient;
  int rounded;
  float new_coord;
  quotient=number/stepsize;
  rounded=(int)quotient;
  new_coord=rounded*stepsize;

  if(min_or_max==1){ // max value
    if(new_coord<=number){
      rounded+=1;
      new_coord=rounded*stepsize;
    }
  }else if(min_or_max==0){ // min value
    if(new_coord>=number){
      rounded-=1;
      new_coord=rounded*stepsize;
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


void get_enrg(map<int,vertex>& grid, vector<atom>& prot, vector<int>& vrtxList){
  int j;
  cout<<endl<< "Searching for potential interactions at each grid intersection"<< endl;
  int cvrtx=0;

  //OLD #pragma omp parallel for shared(prot) schedule(static) num_threads(4)
  #pragma omp parallel for schedule(dynamic) num_threads(8)
  for(int i=0; i<vrtxList.size(); i++){
    try{
      vertex& m=grid.at(vrtxList.at(i));
      cvrtx++;
      
      int flag=0;
      if(m.grid[0]!=1 && m.grid[1]!=1 && m.grid[2]!=1 && m.grid[3]==1){
        continue;
      }
    
      if(printDetails==1){ cout<<endl<<"Vertex id: "<< m.id <<" "<<m.x<<" "<<m.y<<" "<<m.z<<endl; }

      for(int probe=0; probe<nbOfProbes; probe++){ //Iterate each probe
        float enrg_sum=0.00;
        int countAtms=0;

        if(printDetails==1){ cout<<endl<<"### PROBE "<<probe<<" ###"<<endl; }

        for(int j=0; j<prot.size(); j++){ //Iterate each atom for this probe at this grid intersection
          enrg_sum+=calcNrg_FAINT(m,prot.at(j),probe,countAtms);
          printDetails=0;
        }
        if(countAtms>0){
          if(enrg_sum<epsilons[probe] || (fabs(enrg_sum-epsilons[probe]))<0.0001){
            m.ints[probe]=1;
            flag=1;
          }
        }
      }
      if(flag==1){
        for(int gi=0; gi<4; gi++){
          if(m.grid[gi]==1){
            ss[gi]++;
          }
        }
      }
    }catch (const out_of_range& oor) {
      //cerr << "Out of Range error: " << oor.what() << '\n';
    }catch(exception& e){
      cout <<"error"<<endl;
    }
  }
  cout<<"count vrtx in get_enrg: "<<cvrtx<<endl;
  //cout<<"Calculated non-bonded energies at "<<countVrtx<<" grid intersections using "<< nbOfProbes <<" probe(s)."<<endl;
}

// int Grid::smooth(){
//   float xi,yi,zi;
//   float** tmp_sm;
//   map<int,vertex>::iterator it;
//   map<int,vertex>::iterator isit;
//   int probe,i,m,n;
//   float sum_100,sum_50,sum_25;
//   float count_100,count_50,count_25;
//   int count_passed=0;
//   int count_found=0;
//   float smoothing_dist;

//   cout<<endl<<"Smoothing energies"<<endl;
//   for(it=GRID.begin();it!=GRID.end();it++){
//     //printf("\n\nsmoothing: %5.2f, %5.2f, %5.2f",it->second.x,it->second.y,it->second.z);

//     //Allocate temporary memory
//     tmp_sm = new float*[nbOfProbes];
//     for(m=0; m<nbOfProbes; m++){
//       tmp_sm[m] = new float[6];
//     }

//     for(m=0; m<nbOfProbes; m++){
//       for(n=0; n<=5; n++){ //Set em all to 0
//         tmp_sm[m][n]=0.0;
//       }
//     }

//     count_passed=0;
//     count_found=0;

//     //Check in which grid resolution this vertex is in
//     if(it->second.grid[2]==1){
//       smoothing_dist=1.0;

//       //cout<<endl<<"In grid 2.0 "<<it->second.x<< " "<<it->second.y<< " "<< it->second.z<<endl<<"Here are its neibrs:"<<endl;

//       //Go through all the coords in the 1.0A edge cube around this coord
//       for(zi=(it->second.z-smoothing_dist); zi<=(it->second.z+smoothing_dist); zi+=stepsize){
//         for(yi=(it->second.y-smoothing_dist); yi<=(it->second.y+smoothing_dist); yi+=stepsize){
//           for(xi=(it->second.x-smoothing_dist); xi<=(it->second.x+smoothing_dist); xi+=stepsize){
//             count_passed++;

//             i=generateID(width,height,(int)((xi-min_x)/stepsize)+1,(int)((yi-min_y)/stepsize)+1,(int)((zi-min_z)/stepsize)+1);
//             isit = GRID.find(i);

//             //printf("\n\t %6.2f, %6.2f, %6.2f",xi,yi,zi);
//             if(isit != GRID.end()){
//               count_found++;
//               if(stepsize < 0.5){
//                 if(is_coord_in_cube(xi, yi, zi, it->second.x, it->second.y, it->second.z, 0.25)==1){
//                   //printf(" cube 0.25 s<0.5");
//                   for(probe=0; probe<nbOfProbes; probe++){
//                     tmp_sm[probe][0]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][1]+=1.0;
//                     tmp_sm[probe][2]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][3]+=1.0;
//                     tmp_sm[probe][4]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][5]+=1.0;
//                     //printf("\n\t\t[%d]{0.25} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][0],tmp_sm[probe][1]);
//                     //printf("\n\t\t[%d]{0.5} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][2],tmp_sm[probe][3]);
//                     //printf("\n\t\t[%d]{1.0} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][4],tmp_sm[probe][5]);
//                   }
//                 }else if(is_coord_in_cube(xi, yi, zi, it->second.x, it->second.y, it->second.z, 0.5)==1){
//                   //printf(" cube 0.5 s<0.5");
//                   for(probe=0; probe<nbOfProbes; probe++){
//                     tmp_sm[probe][2]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][3]+=1.0;
//                     tmp_sm[probe][4]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][5]+=1.0;
//                     //printf("\n\t\t[%d]{0.5} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][2],tmp_sm[probe][3]);
//                     //printf("\n\t\t[%d]{1.0} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][4],tmp_sm[probe][5]);
//                   }
//                 }else{
//                   //printf(" cube 1.0 s<0.5");
//                   for(probe=0; probe<nbOfProbes; probe++){
//                     tmp_sm[probe][4]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][5]+=1.0;
//                     //printf("\n\t\t[%d]{1.0} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][4],tmp_sm[probe][5]);
//                   }
//                 }
//               }else if(stepsize < 1.0){
//                 if(is_coord_in_cube(xi, yi, zi, it->second.x, it->second.y, it->second.z, 0.5)==1){
//                   //printf(" cube 0.5 s<1.0");
//                   for(probe=0; probe<nbOfProbes; probe++){
//                     tmp_sm[probe][2]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][3]+=1.0;
//                     tmp_sm[probe][4]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][5]+=1.0;
//                     //printf("\n\t\t[%d]{0.5} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][2],tmp_sm[probe][3]);
//                     //printf("\n\t\t[%d]{1.0} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][4],tmp_sm[probe][5]);
//                   }
//                 }else{
//                   //printf(" cube 1.0 s<1.0");
//                   for(probe=0; probe<nbOfProbes; probe++){
//                     tmp_sm[probe][4]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][5]+=1.0;
//                     //printf("\n\t\t[%d]{1.0} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][4],tmp_sm[probe][5]);
//                   }
//                 }
//               }else if(stepsize < 2.0){
//                 //printf(" cube 1.0 s<2.0");
//                 for(probe=0; probe<nbOfProbes; probe++){
//                   tmp_sm[probe][4]+=isit->second.nrg[probe][0];
//                   tmp_sm[probe][5]+=1.0;
//                   //printf("\n\t\t[%d]{1.0} %8.6f sum: %8.6f count: %8.6f",probe,isit->second.nrg[probe][0],tmp_sm[probe][4],tmp_sm[probe][5]);
//                 }
//               }
//             }else{
//               //cout<<" coord does not exist"<<endl;
//             }
//           }
//         }
//       }
//     }else if(it->second.grid[1]==1 && stepsize < 1.0){
//       smoothing_dist=0.5;
//       //cout<<endl<<"In grid 1.0 "<<it->second.x<< " "<<it->second.y<< " "<< it->second.z<<endl<<"Here are its neibrs:"<<endl;

//       //Go through all the coords in the 0.5A edge cube around this coord
//       for(zi=(it->second.z-smoothing_dist); zi<=(it->second.z+smoothing_dist); zi+=stepsize){
//         for(yi=(it->second.y-smoothing_dist); yi<=(it->second.y+smoothing_dist); yi+=stepsize){
//           for(xi=(it->second.x-smoothing_dist); xi<=(it->second.x+smoothing_dist); xi+=stepsize){
//             count_passed++;

//             i=generateID(width,height,(int)((xi-min_x)/stepsize)+1,(int)((yi-min_y)/stepsize)+1,(int)((zi-min_z)/stepsize)+1);
//             isit = GRID.find(i);
            
//             if(isit != GRID.end()){
//               count_found++;

//               if(stepsize < 0.5){
//                 if(is_coord_in_cube(xi, yi, zi, it->second.x, it->second.y, it->second.z, 0.25)==1){
//                   //printf(" cube 0.25");
//                   for(probe=0; probe<nbOfProbes; probe++){
//                     tmp_sm[probe][0]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][1]+=1.0;
//                     tmp_sm[probe][2]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][3]+=1.0;
//                     //printf("\n\t\t[%d]{0.5} %8.6f sum: %8.6f count: %8.6f",probe,isit->second->nrg[probe][0],tmp_sm[probe][2],tmp_sm[probe][3]);
//                     //printf("\n\t\t[%d]{0.25} %8.6f sum: %8.6f count: %8.6f",probe,isit->second->nrg[probe][0],tmp_sm[probe][0],tmp_sm[probe][1]);
//                   }
//                 }else{
//                   //printf(" cube 0.5");
//                   for(probe=0; probe<nbOfProbes; probe++){
//                     tmp_sm[probe][2]+=isit->second.nrg[probe][0];
//                     tmp_sm[probe][3]+=1.0;
//                     //printf("\n\t\t[%d]{0.5} %8.6f sum: %8.6f count: %8.6f",probe,isit->second->nrg[probe][0],tmp_sm[probe][2],tmp_sm[probe][3]);
//                   }
//                 }
//               }else {
//                 //printf(" cube 0.5");
//                 for(probe=0; probe<nbOfProbes; probe++){
//                   tmp_sm[probe][2]+=isit->second.nrg[probe][0];
//                   tmp_sm[probe][3]+=1.0;
//                   //printf("\n\t\t[%d]{0.5} %8.6f sum: %8.6f count: %8.6f",probe,isit->second->nrg[probe][0],tmp_sm[probe][2],tmp_sm[probe][3]);
//                 }
//               }
//             }else{
//               //cout<<"coord does not exist"<<endl;
//             }
//           }
//         }
//       }
//     }else if(it->second.grid[0]==1 && stepsize<0.5){ //If it's in 0.5 we'll need to smooth for 0.25A
//       smoothing_dist=0.25;
//       //printf(" [0.5]\nsmooth_dist %7.3f\n\tNeibrs:",smoothing_dist);

//       //Go through all the coords in the 0.25A edge cube around this coord
//       for(zi=(it->second.z-smoothing_dist); zi<=(it->second.z+smoothing_dist); zi+=stepsize){
//         for(yi=(it->second.y-smoothing_dist); yi<=(it->second.y+smoothing_dist); yi+=stepsize){
//           for(xi=(it->second.x-smoothing_dist); xi<=(it->second.x+smoothing_dist); xi+=stepsize){
//             count_passed++;

//             i=generateID(width,height,(int)((xi-min_x)/stepsize)+1,(int)((yi-min_y)/stepsize)+1,(int)((zi-min_z)/stepsize)+1);
//             isit = GRID.find(i);
            
//             //printf("\n\t %6.2f, %6.2f, %6.2f",xi,yi,zi);
//             if(isit != GRID.end()){
//               count_found++;
//               //printf(" cube 0.25");
//               for(probe=0; probe<nbOfProbes; probe++){
//                 tmp_sm[probe][0]+=isit->second.nrg[probe][0];
//                 tmp_sm[probe][1]+=1.0;
//                 //printf("\n\t\t[%d]{0.25} %8.6f sum: %8.6f count: %8.6f",probe,isit->second->nrg[probe][0],tmp_sm[probe][0],tmp_sm[probe][1]);
//               }
//             }else{
//               //cout<<"coord does not exist"<<endl;
//             }
//           }
//         }
//       }
//       //cout<<"passed "<<count_passed<<endl;
//     }else{

//     }

//     //Calculate and store the smoothed probes
//     //printf("\nPASSED %d FOUND %d",count_passed,count_found);
//     for(probe=0; probe<nbOfProbes; probe++){
//       if(tmp_sm[probe][1]>0.00){
//         it->second.nrg[probe][1]=tmp_sm[probe][0]/tmp_sm[probe][1];
//         it->second.pb_bool[probe][1]=1;
//       }else{
//         it->second.pb_bool[probe][1]=0; //No value = well assign 0
//       }
//       if(tmp_sm[probe][3]>0.00){
//         it->second.nrg[probe][2]=tmp_sm[probe][2]/tmp_sm[probe][3];
//         it->second.pb_bool[probe][2]=1;
//       }else{
//         it->second.pb_bool[probe][2]=0; //No value = well assign 0
//       }
//       if(tmp_sm[probe][5]>0.00){
//         it->second.nrg[probe][3]=tmp_sm[probe][4]/tmp_sm[probe][5];
//         it->second.pb_bool[probe][3]=1;
//       }else{
//         it->second.pb_bool[probe][3]=0; //No value = well assign 0
//       }
//       //printf("\n[0] %15.10f [1] %15.10f [2] %15.10f [3] %15.10f [4] %15.10f [5] %15.10f",tmp_sm[probe][0],tmp_sm[probe][1],tmp_sm[probe][2],tmp_sm[probe][3],tmp_sm[probe][4],tmp_sm[probe][5]);
//       //printf("\nFINAL SMOOTHED ENERGY [probe: %d] [0.25] %22.15f [0.5] %22.15f [1.0] %22.15f\n",probe,it->second.nrg[probe][1],it->second.nrg[probe][2],it->second.nrg[probe][3]);
//     }
    
//     //printf("\nDeleting memory\n");
//     for (m=0; m<nbOfProbes; m++){
//       delete [] tmp_sm[m];
//     }
//     delete [] tmp_sm;  
//   }
//   return(0);
// }

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

void Grid::writeMif(vector<atom>& prot){
  map<int,vertex>::iterator it;

  FILE* fpNew;
  int probe;
  float d=0.0;
  int bsFlag=0;
  int i,j;
  string na="NA";
  string mifFileNew=outBase + tag + ".mif";

  cout<<endl<< "Writing Mif File"<<endl;

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
  fprintf(fpNew,"#nbOfProbes %2d\n",nbOfProbes);
  for(i=0; i<nbOfProbes; i++){
    fprintf(fpNew,"#probe[%d] %4d -> %3s %3s %5s\n",i,probesList[i],atomSingle->res[probesList[i]].c_str(),atomSingle->atom[probesList[i]].c_str(),atSingle->name[atomSingle->atId[probesList[i]]].c_str());
  }
  fprintf(fpNew,"#stepsize %4.2f\n",stepsize);
  fprintf(fpNew,"#atom_probe_distance_threshold %7.4f\n",atmPbMaxDist);
  fprintf(fpNew,"#protein_grid_distance %7.4f to %7.4f\n",sqrt(minGridDist),sqrt(maxGridDist));
  fprintf(fpNew,"#grid_width %d\n",width);
  fprintf(fpNew,"#grid_height %d\n",height);
  fprintf(fpNew,"#gs %d %d %d %d\n",vrtx200,vrtx150,vrtx100,vrtx050);
  fprintf(fpNew,"#ss %d %d %d %d\n",ss[0],ss[1],ss[2],ss[3]);

  int cvrtx=0;
  for(it=GRID.begin();it!=GRID.end();it++){
    fprintf(fpNew, "%7.2f %7.2f %7.2f",
      it->second.x,
      it->second.y,
      it->second.z);
    for(probe=0; probe<nbOfProbes; probe++){
      fprintf(fpNew, " %1d",it->second.ints[probe]);
    }
    for(i=0; i<4; i++){
      fprintf(fpNew, " %1d",it->second.grid[i]);
    }
    fprintf(fpNew, "\n");
    cvrtx++;
  }

  if(printAtoms==1){
    for(i=0; i<prot.size(); i++){
      if(prot[i].h!=1) continue;
      if(prot[i].bs==1) continue;
      bsFlag=0;
      for(it=GRID.begin();it!=GRID.end();it++){
        d=dist_3d(prot[i].x,prot[i].y,prot[i].z,it->second.x,it->second.y,it->second.z);
        if(d<caT || fabs(d-caT)<0.001){
          bsFlag=1;
          break;
        }
      }
      if(bsFlag==1){
        // cout<<endl<<prot[i].resn<<" "<<prot[i].resnb<<" "<<prot[i].atomn<<" "<<prot[i].chain;
        for(j=0; j<prot.size(); j++){
            if(prot[j].resn.compare(prot[i].resn)==0 && prot[j].resnb==prot[i].resnb && prot[j].chain.compare(prot[i].chain)==0){
              // cout<<" found "<<prot[j].resn<<" "<<prot[j].resnb<<" "<<prot[j].atomn<<" "<<prot[j].chain;
              prot[j].bs=1;
            }
        }
      }
    }  
    for(j=0; j<prot.size(); j++){
      fprintf(fpNew,"#ATOM %3s %4d %4s %5d %s %8.3f %8.3f %8.3f %d %d\n",prot[j].resn.c_str(),prot[j].resnb,prot[j].atomn.c_str(),prot[j].atomnb,prot[j].chain.c_str(),prot[j].x,prot[j].y,prot[j].z,prot[j].mif,prot[j].bs);
    }   
  }


  fclose(fpNew);
  cout<<"count vrtx file "<<cvrtx<<endl;
}

void stripSpace(string &str) {
  for (int i=0;i<str.length();i++) 
    if (str[i]==' ') {
    str.erase(i,1);
    i--;
  }
}

void storeAtomSingle_FAINT(atomSingleStruc* atomSingle){
  string file= basePath + "/forcefield_files/atomSingle.isomif";
  ifstream infile(file.c_str());
  string line;
  string res, atom;
  int atomId,atId;

  while(getline(infile,line)){
    if(line.compare(0,1,"#")!=0){ nbOfAtoms++; }
  }

  infile.clear();
  infile.seekg(0, ios::beg);
  atomSingle->count=nbOfAtoms;
  atomSingle->atId=new int[nbOfAtoms];
  atomSingle->res=new string[nbOfAtoms];
  atomSingle->atom=new string[nbOfAtoms];
  //cout<< "Getting AtomSingle"<<endl;
  while(getline(infile,line)){
    if(line.compare(0,1,"#")!=0){
      std::stringstream test(line);
      test >> atomId >> res >> atom >> atId;
      atomSingle->atId[atomId]=atId;
      atomSingle->res[atomId]=res;
      atomSingle->atom[atomId]=atom;
      if(printDetails==1){ cout<< atomId<< " "<< res<< " "<< atom<< " "<< atId<< endl; }
    }
  }
}

void storeAtSingle_FAINT(atSingleStruc* atSingle){
  string file= basePath + "/forcefield_files/atSingle.isomif";
  ifstream infile(file.c_str());
  string line, name, pseudo;
  int HbD, HbA, aromatic, charged, atId, HydroPhb;

  while(getline(infile,line)){
    if(line.compare(0,1,"#")!=0){ nbOfAts++; }
  }

  atSingle->name=new string[nbOfAts];
  atSingle->HbD_FAINT=new int[nbOfAts];
  atSingle->HbA_FAINT=new int[nbOfAts];
  atSingle->aromatic_FAINT=new int[nbOfAts];
  atSingle->charged_FAINT=new int[nbOfAts];
  atSingle->hydrophobic_FAINT=new int[nbOfAts];

  infile.clear();
  infile.seekg(0, ios::beg);
  //printf("%42s %5s %5s %5s %5s\n","Probe Name","HbD","HbA","Arom","Chrgd");
  while(getline(infile,line)){
    if(line.compare(0,1,"#")!=0){
      std::stringstream test(line);
      test >> atId >> name >> pseudo >> HbD >> HbA >> aromatic >> charged >> HydroPhb;
      atSingle->name[atId]=name + " " + pseudo;
      atSingle->HbD_FAINT[atId]=HbD;
      atSingle->HbA_FAINT[atId]=HbA;
      atSingle->aromatic_FAINT[atId]=aromatic;
      atSingle->charged_FAINT[atId]=charged;
      atSingle->hydrophobic_FAINT[atId]=HydroPhb;
      if(printDetails==1){ 
        cout<< atId<< " "<< name<< " "<< pseudo << " "<< HbD<< " "<< HbA<< " "<< aromatic<< " "<< charged<<" "<<HydroPhb<< endl;
      }
    }
  }
}

void storeAtPair_FAINT(atPairStruc* atPair){
  string line;
  int id1,id2;
  float epsilonFAINT;

  atPair->epsilon_FAINT=new float[nbOfAts*nbOfAts];
  string atPairFile= basePath + "/forcefield_files/atPair.isomif";
  ifstream infile(atPairFile.c_str());
  while(getline(infile,line)){
    std::stringstream test(line);
    test >> id1 >> id2 >> epsilonFAINT;
    atPair->epsilon_FAINT[(id1*nbOfAts)+id2]=epsilonFAINT;
    atPair->epsilon_FAINT[(id2*nbOfAts)+id1]=epsilonFAINT;
    if(printDetails==1){ cout<< id1 << " "<< id2 << " "<< atPair->epsilon_FAINT[(id1*nbOfAts)+id2]<< endl; }
  }
}

int getProbes_FAINT(){
  string file= basePath + "/forcefield_files/probes.isomif";
  ifstream infile(file.c_str());
  string line;
  int count=0;
  int id;
  float distTmin;
  float distTmax;
  float epsilon;
  while(getline(infile,line)){
    nbOfProbes++;
  }
  infile.clear();
  infile.seekg(0, ios::beg);
  probesList = new int[nbOfProbes];
  pbDistTmin = new float[nbOfProbes];
  pbDistTmax = new float[nbOfProbes];
  epsilons = new float[nbOfProbes];

  if(printDetails==1){ cout<< endl<< "Getting probes"<< endl; }
  while(getline(infile,line)){
    std::stringstream test(line);
    test >> id >> distTmin >> distTmax >> epsilon;
    probesList[count]=id;
    pbDistTmin[count]=distTmin;
    pbDistTmax[count]=distTmax;
    epsilons[count]=epsilon;
    // cout<<count<<" "<<epsilon<<endl;
    if(printDetails==1){ cout<< id<< " " <<atomSingle->atId[id]<< " "<< atomSingle->atom[id]<< " "<< atomSingle->res[id]<< endl; }
    count++;
  }
  if(printDetails==1){ cout << "Nb of probes: " << nbOfProbes << endl; }
  return(nbOfProbes);
}

double calcNrg_FAINT(vertex& vrtx, atom& atm, int pbId, int& count_atoms){
  float dist,alpha,epsilon,angle;
  int atomAtId,pbAtId;
  double energy;
  energy=0.0;
  alpha=1.0;
  angle=1.0;
  float angleThresh=40.00;
  int tVrtxId=-1;
  float rDist,rpDist;

  if(atm.atomn.compare("N")==0 && pbId==3 && 0){
    printDetails=1;
  }else{
    printDetails=0;
  }

  if(atm.resnb==518 && atm.atomn.compare("CZ")==0 && pbId==5){
    printDetails=1;
  }

  //Get distance between probe and atom
  dist=dist_3d(atm.x,atm.y,atm.z,vrtx.x,vrtx.y,vrtx.z);
  if(dist > pbDistTmax[pbId] || dist < pbDistTmin[pbId] || dist > atmPbMaxDist){
    if(printDetails==1){
      // cout<<"Dist "<<dist<<" - Atom too far or too close to probe."<<endl;
    }
    return(energy);
  }else{

    atomAtId=atomSingle->atId[atm.atomId];
    pbAtId=atomSingle->atId[probesList[pbId]];
    epsilon=atPair->epsilon_FAINT[(atomAtId*nbOfAts)+pbAtId];

    if(printDetails==1){
      cout<< endl <<atm.resn<< " "<<atm.resnb<< " "<< atm.atomn<<" "<< atomAtId<< " "<< atSingle->name[atomAtId]<< " x "<< atm.x<<" y "<<atm.y<<" z "<< atm.z<<" xr" << atm.xr<<" yr "<<atm.yr<<" xr "<< atm.zr<<endl;
      cout<< "Probe "<< pbId<<" "<< atSingle->name[pbAtId] <<" "<< pbAtId <<" id: "<<vrtx.id<<" "<< vrtx.x<<" "<<vrtx.y<<" "<<vrtx.z<<endl;
    }

    //Hbond Donnor/Acceptor
    if((atSingle->HbD_FAINT[atomAtId]==1 && atSingle->HbA_FAINT[pbAtId]==1) || (atSingle->HbA_FAINT[atomAtId]==1 && atSingle->HbD_FAINT[pbAtId]==1)){
      alpha=1.0;

      rDist=dist_3d(atm.x,atm.y,atm.z,atm.xr,atm.yr,atm.zr);
      rpDist=dist_3d(vrtx.x,vrtx.y,vrtx.z,atm.xr,atm.yr,atm.zr);

      if(printDetails==1){
        cout <<"Hbond Acceptor/Donor "<<endl;
        cout<< "Dist "<< dist << " rDist "<< rDist<< " rpDist "<< rpDist<<" rDir "<< atm.rDir <<endl;
      }

      angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
      angle=acos(angle)* 180.0 / PI;
      if(atm.rDir==0){ angle=180.00-angle; }

      if(atm.resn.compare("ASN")==0 && atm.atomn.compare("OD1")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("ASN")==0 && atm.atomn.compare("ND2")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("GLN")==0 && atm.atomn.compare("OE1")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("GLN")==0 && atm.atomn.compare("NE2")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("SER")==0 && atm.atomn.compare("OG")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("THR")==0 && atm.atomn.compare("OG1")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("ARG")==0 && atm.atomn.compare("NH1")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("ARG")==0 && atm.atomn.compare("NH2")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("LYS")==0 && atm.atomn.compare("NZ")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("ASP")==0 && atm.atomn.compare("OD1")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("ASP")==0 && atm.atomn.compare("OD2")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("GLU")==0 && atm.atomn.compare("OE1")==0){
        angleThresh=90.00;
      }else if(atm.resn.compare("GLU")==0 && atm.atomn.compare("OE2")==0){
        angleThresh=90.00;
      }

      if(printDetails==1){
        cout <<"Angle "<< angle << " thresh " << angleThresh << " | ref atom: "<<atm.xr<<" "<<atm.yr<<" "<<atm.zr << " dir: "<<atm.rDir<<endl;
      }

      if(angle>angleThresh){
        if(printDetails==1){ cout<<"Angle over threshold"<<endl; }
        return(energy);
      }
    }else if(atSingle->aromatic_FAINT[atomAtId]==1 && atSingle->aromatic_FAINT[pbAtId]==1){ //aromatic interaction
      alpha=0.75;
      rDist=dist_3d(atm.x,atm.y,atm.z,atm.xr,atm.yr,atm.zr);
      rpDist=dist_3d(vrtx.x,vrtx.y,vrtx.z,atm.xr,atm.yr,atm.zr);
      angle=(pow(dist,2.0)+pow(rDist,2.0)-pow(rpDist,2.0))/(2*dist*rDist);
      angle=acos(angle)* 180.0 / PI;
      if(angle>90 || fabs(180.00-angle) < 0.001){
        angle=180-angle;
      }

      if(printDetails==1){
        cout<< "aromatic"<<endl;
        cout<< "Angle "<< angle<<endl;
      }

      if(angle>40.00 && angle < 80.00){
        if(printDetails==1){ cout <<"angle between 60 and 80"<<endl; }
        return(energy);
      }else{
        if(printDetails==1){ cout <<"angle < 60 or > 80"<<endl; }
      }

    }else if(atSingle->charged_FAINT[atomAtId]==1 && atSingle->charged_FAINT[pbAtId]==1){ //charged interaction
      if(printDetails==1){ cout << "charged couple"<< endl; }
      alpha=0.75;
    }else if(atSingle->hydrophobic_FAINT[pbAtId]==1){//if its a hydrophoic probe
      if(printDetails==1){ cout << "Hydrophobic probe"<< endl; }
      alpha=0.75;
    }

    count_atoms++;
    energy=(epsilon)*(exp(-1.0*dist*alpha));

    if(printDetails==1){
      cout<< "epsilon: " << epsilon<< " alpha: "<< alpha<< " dist: "<< dist<< " -> NRG: "<< energy<< endl<<endl;
    }
    return(energy);
  }
}