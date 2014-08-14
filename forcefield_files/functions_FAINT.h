int getProbes_FAINT();
void storeAtPair_FAINT(atPairStruc*);
void storeAtomSingle_FAINT(atomSingleStruc*);
void storeAtSingle_FAINT(atSingleStruc*);
double calcNrg_FAINT(vertex&, atom&, int, int&);

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