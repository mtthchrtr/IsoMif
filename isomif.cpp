#include "isomif.h"

/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int main(int argc, char *argv[]){
  char tmp[2000];
  FILE* fp;
  int i,j,k,l,cg,cgs;
  int count=0;
  vector<node>::iterator vit1;
  pNode agv=NULL;
  // float cen_a[3];
  // float cen_b[3];
  // gsl_matrix *mat_r  = gsl_matrix_alloc(3,3);

  if(read_commandline(argc, argv)==24){
    printf("\nread cmd line not OK\n");
    return(24);
  }

  setJTT(jttt);

  if(pairwiseF.compare("")!=0){
    getPairwise();
    open_file_ptr(&fpout,out_file,1);
  }else{
    pwRun npw;
    npw.mif1= nrg_file1;
    npw.mif2= nrg_file2;
    npw.rnc1=rnc1;
    npw.rnc2=rnc2;
    pw.push_back(npw);
  }

  for(int pwr=0; pwr<pw.size(); pwr++){
  
  cout<<pwr<<" mif1 "<<pw[pwr].mif1<<" mif2 "<<pw[pwr].mif2<<" rnc1 "<<pw[pwr].rnc1<<" rnc2 "<<pw[pwr].rnc2<<endl;
  nrg_file1=pw[pwr].mif1;
  nrg_file2=pw[pwr].mif2;
  rnc1=pw[pwr].rnc1;
  rnc2=pw[pwr].rnc2;
  topT=-1.0;
  topN=-1;

  if(get_info(nrg_file1,nrg_file2)==24){ return(24); }

  if(emptOut==1){
    sprintf(tmp,"REMARK command: %s\nREMARK wsimfn: %d\nREMARK tag1: %s\nREMARK tag2: %s\nREMARK rnc1: %s\nREMARK rnc2: %s\n",cmdLine,wrfn,tag1.c_str(),tag2.c_str(),rnc1.c_str(),rnc2.c_str());  
  }else{
    sprintf(tmp,"REMARK command: %s\nREMARK mif_file_1: %s\nREMARK mif_file_2: %s\nREMARK nb_of_probes: %d\nREMARK C-alpha_dDist: %5.2f\nREMARK pseudocenter_dDist: %5.2f\nREMARK dDist: %5.2f\nREMARK jtt_threshold: %d\nREMARK max_nodes: %d\nREMARK commont int : %d\nREMARK wsimfn: %d\nREMARK tag1: %s\nREMARK tag2: %s\nREMARK rnc1: %s\nREMARK rnc2: %s\n",cmdLine,nrg_file1.c_str(),nrg_file2.c_str(),nb_of_probes,ca_dDist,ps_dDist,dDist,jttt,maxNodes,commonInt,wrfn,tag1.c_str(),tag2.c_str(),rnc1.c_str(),rnc2.c_str());  
  }
  
  strcpy(outH,tmp);

  //Print info to User
  // printf("--# Get Info #--\n\nEnergy file 1: %s\nEnergy file 2: %s\nOutfile: %s\nNumber of probes: %d\n\n--# Parameter Values #--\n\nC-alpha_dDist: %5.2f\nREMARK pseudocenter_dDist: %5.2f\nNode dDist: %5.2f\nneibr dDist: %5.2f\nJTT threshold: %d\nMaximum nodes: %d\nGet all C-alpha cliques: %d\nMinimum number of Common interaction: %d\nPrint details: %d\n",nrg_file1.c_str(),nrg_file2.c_str(),out_file,nb_of_probes,ca_dDist,ps_dDist,dDist,neibr_dDist,jttt,maxNodes,bkAll,commonInt,printDetails);
  // cout<<"CG begin "<<steps.front()<<" cg_start "<<cg_start<<" end "<<steps.back()<<endl<<endl;
  // cout<<rnc1<<" "<<rnc2<<endl;
  createVrtxVec(nrg_file1,mif1,prot1,ss1,ss1m,caSize1,pseudoL1,rnc1,lig1);
  createVrtxVec(nrg_file2,mif2,prot2,ss2,ss2m,caSize2,pseudoL2,rnc2,lig2);
  cout<<"mif1 size: "<<mif1.size()<<endl;
  cout<<"mif2 size: "<<mif2.size()<<endl;

  // for(int i=0; i<lig1.size(); i++){
  //   cout<<lig1[i].resn<<" "<<lig1[i].resnb<<" "<<lig1[i].chain<<" "<<lig1[i].atomn<<" "<<lig1[i].coor[0]<<" "<<lig1[i].coor[1]<<" "<<lig1[i].coor[2]<<endl;
  // }

  // for(int i=0; i<lig2.size(); i++){
  //   cout<<lig2[i].resn<<" "<<lig2[i].resnb<<" "<<lig2[i].chain<<" "<<lig2[i].atomn<<" "<<lig2[i].coor[0]<<" "<<lig2[i].coor[1]<<" "<<lig2[i].coor[2]<<endl;
  // }
  // return(0);

  //Set which vertices to be considered if cg_start > 0
  if(cg_start>-1){
    for(i=0; i<mif1.size(); ++i){
      mif1.at(i).cg[cg_start]=1;
    }
    for(i=0; i<mif2.size(); ++i){
      mif2.at(i).cg[cg_start]=1;
    }
  }

  cout << "--# Starting coarsegrain steps #--\n";

  //Start the coarsegrain steps
  for(int cs=0; cs<steps.size(); cs++){
    nCliques=0;
    cout<<endl<<"Coarse-Grain Step "<<steps[cs]<<endl;

    //Step 2, superimpose mifs using rotation matrix derived from atom list superimposition
    if(steps[cs]==-2){
      vector<float> la;
      vector<float> lb;

      if(list1.size()==0){
        cout<<"You must provide two lists of corresponding atom IDs to superimpose for stage -2 using argument -q. Ex: -q 0,1,2,3 38,46,47,53"<<endl;
      }

      //Create coords lists
      for(int i=0; i<list1.size(); i++){
        for(int j=0; j<prot1.size(); j++){
          if(prot1[j].atomnb==list1[i]){
            for(int c=0; c<3; c++){
              la.push_back(prot1[j].coor[c]);
            }
            break;
          }
        }
        for(int j=0; j<prot2.size(); j++){
          if(prot2[j].atomnb==list2[i]){
            for(int c=0; c<3; c++){
              lb.push_back(prot2[j].coor[c]);
            }
            break;
          }
        }
      }

      if(la.size()!=lb.size()) cout<<"The -2 superimposition requires to have same number of coordinates."<<endl;

      Clique nc;
      nc.cg=-2;
      nc.mat_r=gsl_matrix_alloc(3,3);
      for(int i=0; i<3; i++) {
        nc.cen_a[i]=0.0;
        nc.cen_b[i]=0.0;
      }
      gsl_matrix_set_zero(nc.mat_r);

      nc.det=calcRot(la,lb,nc.cen_a,nc.cen_b,nc.mat_r);

      // for(int i=0; i<3; i++) {
      //   cout<<"cen_a "<<i<<" "<<nc.cen_a[i]<<endl;
      //   cout<<"cen_b "<<i<<" "<<nc.cen_b[i]<<endl;
      // }

      // for(int i=0; i<3; i++) {
      //   for(int j=0; j<3; j++) {
      //     cout<<i<<" "<<j<<" "<<gsl_matrix_get(nc.mat_r,i,j)<<endl;
      //   }
      // }

      // cout<<"Rotating vertexes of Mif 1 onto Mif 2 using list of atoms..."<<endl;
      //Rotate mif 1 onto mif 2
      for(int v=0; v<mif1.size(); v++){
        for(int i=0; i<3; i++){
          mif1[v].ncoor[i]=nc.cen_b[i];
          for(int j=0; j<3; j++){
            mif1[v].ncoor[i]+=(mif1[v].coor[j]-nc.cen_a[j])*gsl_matrix_get(nc.mat_r,i,j);
          }
        }
      }

      // cout<<"Finding corresponding vertexes..."<<endl;
      float dist=0.0;
      for(int u=0; u<mif1.size(); u++){
        if(mif1[u].grid[cg2]!=1) continue;
        for(int v=0; v<mif2.size(); v++){
          if(mif2[v].grid[cg2]!=1) continue;
          dist=dist3d(mif1[u].ncoor,mif2[v].coor);        
          if(dist < dDist || fabs(dist-dDist)<0.001){ //If passes distance threshold
            for(int i=0; i<nb_of_probes; i++){
              if(mif1[u].pb[i]==1 && mif2[v].pb[i]==1){
                // cout<<i<<" - "<<mif1[u].ncoor[0]<<" "<<mif1[u].ncoor[1]<<" "<<mif1[u].ncoor[2]<<" "<<mif2[v].coor[0]<<" "<<mif2[v].coor[1]<<" "<<mif2[v].coor[2]<<" - "<<mif1[u].pb[i]<<" "<<mif2[v].pb[i]<<endl;  
                mif1[u].m[i]=1;
                mif2[v].m[i]=1;
              }
            }
          }
        }
      }

      for(int u=0; u<mif1.size(); u++){
        for(int i=0; i<nb_of_probes; i++){
          if(mif1[u].m[i]==1){
            nc.va.push_back(mif1[u]);
            break;
          }
        }
      }
      for(int v=0; v<mif2.size(); v++){
        for(int i=0; i<nb_of_probes; i++){
          if(mif2[v].m[i]==1){
            nc.vb.push_back(mif2[v]);
            break;
          }
        }
      }

      nc.nbNodes=nc.va.size()+nc.vb.size();
      nc.tani=( ((float)nc.va.size()/(float)ss1[cg2]) + ((float)nc.vb.size()/(float)ss2[cg2]) ) / 2.0;
      // cout<<"NEW CLIQUE CG "<<steps[cs]<<" NODES "<<nc.nbNodes<<" TANI "<<nc.tani<<endl;
      cliques.push_back(nc);
    }else{
      bool* conn=NULL;
      vector<node> graph;

      if(steps[cs]!=steps[cs-1]){
        cout<<"Resetting top clique score"<<endl;
        topT=-1.0;
        topN=-1;
      }
      cg=steps[cs];

      //If its not the first stage and its a different grid resolution than previous stage
      //rotate the vertexes using the previous rotation matrix
      if(cs>0 && steps[cs]!=steps[cs-1]){

        //Print Matrix and centers
        // for(int i=0; i<3; i++) {
        //   cout<<"cen_a "<<i<<" "<<cliques.back().cen_a[i]<<endl;
        //   cout<<"cen_b "<<i<<" "<<cliques.back().cen_b[i]<<endl;
        // }
        // for(int i=0; i<3; i++) {
        //   for(int j=0; j<3; j++) {
        //     cout<<i<<" "<<j<<" "<<gsl_matrix_get(cliques.back().mat_r,i,j)<<endl;
        //   }
        // }

        // cout<<"Rotating Mif 1 onto Mif 2 using previous stage..."<<endl;
        for(int v=0; v<mif1.size(); v++){
          for(int i=0; i<3; i++){
            mif1[v].ncoor[i]=cliques.back().cen_b[i];
            for(int j=0; j<3; j++){
              mif1[v].ncoor[i]+=(mif1[v].coor[j]-cliques.back().cen_a[j])*gsl_matrix_get(cliques.back().mat_r,i,j);
            }
          }
        }
        cgs=1; //Wee need to check distance between rotated vertex and remaining vrtx     
      }else{
        cgs=0;
        // cout<<"Not rotating vertexes"<<endl;
      }

      //Create nodes
      cout<<"Creating nodes..."<<endl;
      createNodes(cg,graph,cgs);
      numNodes=graph.size();
      cout<<"NbNodes "<<numNodes<<endl;
    // return(0);
      // cout<<"g.size: "<<graph.size()<<endl;

      // for(i=0; i<graph.size(); i++){
      //   cout<< graph.at(i).uid << " "<< graph.at(i).id << " " << graph.at(i).similarity <<endl;
      // }

      //Sort nodes by similarity
      // sort(graph.begin(), graph.end(), &compareSim);

      cout << "Graph has "<<graph.size() << " nodes."<< endl;

      //If there is too much nodes, sort the list by similarity and keep the max num of nodes with the best similarity
      if(graph.size()>maxNodes){
        int extra=0;
        extra=graph.size()-maxNodes;
        cout<<"Too much nodes must detele some."<<endl<< "There is "<<extra<<" extra nodes"<<endl;
        for(i=0; i<extra; ++i){ graph.pop_back(); }
        numNodes=maxNodes;
        cout << "Graph was shrinked to "<<graph.size() << " nodes."<< endl;
      }

      //reset node ids from 0 to maxNodes
      // for(i=0; i<graph.size(); ++i){
      //   graph.at(i).id=i;
      //   cout<<j<<" id: "<<graph.at(j).id<<endl;
      // }

      //Create adjacency matrix
      adjMat(graph,conn,cg);
      cout<<"numEdges: "<<numEdges<<endl;

      //Find cliques
      cout<<"Entering Bron Kerbosch"<<endl;
      bk(cg,graph,conn);

      //Print nodes in the output file
      clearStep(cg);

      //Delete graph and adjacency matrix
      delete[] conn;
      conn=NULL;

      for(vit1=graph.begin(); vit1<graph.end(); ++vit1){ graph.erase(vit1); }
      vector<node>().swap(graph);
    }
  }

  printNodes();

  for(int v=0; v<mif1.size(); v++){
    mif1[v].nrg.clear();
    mif1[v].ang.clear();
    mif1[v].pb.clear();
    mif1[v].m.clear();
  }
  for(int v=0; v<mif2.size(); v++){
    mif2[v].nrg.clear();
    mif2[v].ang.clear();
    mif2[v].pb.clear();
    mif2[v].m.clear();
  }

  mif1.clear();
  // vector<vertex>().swap(mif1);
  mif2.clear();
  // vector<vertex>().swap(mif2);
  prot1.clear();
  // vector<atom>().swap(prot1);
  prot2.clear();
  // vector<atom>().swap(prot2);
  caSize1=0;
  caSize2=0;
  pseudoL1.clear();
  // vector<pseudoC>().swap(pseudoL1);
  pseudoL2.clear();
  // vector<pseudoC>().swap(pseudoL2);
  lig1.clear();
  // vector<atom>().swap(lig1);
  lig2.clear();
  // vector<atom>().swap(lig2);
  cliques.clear();
  vector<Clique>().swap(cliques);
  for(int k=0; k<4; k++){ ss1[k]=0; ss2[k]=0; ss1m[k]=0; ss2m[k]=0; }

  cout<< "Finished printing nodes and clearing"<<endl;

  }
  fclose(fpout);
  return(0);
}/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void adjMat(vector<node> &graph, bool* &conn, int cg){
  float dist;
  int i=0;
  int j=0;

  long long vsize = (( ( (long long)numNodes * (long long)numNodes ) - (long long)numNodes ) / (long long)2) + (long long)numNodes;
  conn=new bool[vsize];
  cout<<"Size of Connected matrix: "<<vsize<<endl;

  for(long long z=0; z<vsize; z++){ conn[z]=0; }

  numEdges=0;
  for(i=0; i<graph.size(); i++){
    conn[ConnID(i,i)]=1; //Diagonal == 1
    for(j=0; j<i; j++){
      if(cg==-1){
        if(graph.at(i).ca->id != graph.at(j).ca->id && graph.at(i).cb->id != graph.at(j).cb->id){
          dist=fabs(dist3d(graph.at(i).ca->coor,graph.at(j).ca->coor)-dist3d(graph.at(i).cb->coor,graph.at(j).cb->coor));        
          if(dist < ca_dDist || fabs(dist-ca_dDist)<0.001){ //If passes distance threshold
            // cout<<endl<<graph.at(i).a->coor[0]<<" "<<graph.at(i).a->coor[1]<<" "<<graph.at(i).a->coor[2]<<" "<<graph.at(j).a->coor[0]<<" "<<graph.at(j).a->coor[1]<<" "<<graph.at(j).a->coor[2];
            // cout<<endl<<graph.at(i).b->coor[0]<<" "<<graph.at(i).b->coor[1]<<" "<<graph.at(i).b->coor[2]<<" "<<graph.at(j).b->coor[0]<<" "<<graph.at(j).b->coor[1]<<" "<<graph.at(j).b->coor[2];
            // cout<<endl<<"Dist: "<<dist<<endl<<i<<" "<<j<< " ConnID: "<< ConnID(i,j)<<" Connection!"<<endl;
            conn[ConnID(i,j)]=1;
            numEdges++; //create the new edge
          }
        }
      }else if(cg==-3){
        if(graph.at(i).pa->id != graph.at(j).pa->id && graph.at(i).pb->id != graph.at(j).pb->id){
          dist=fabs(dist3d(graph.at(i).pa->coor,graph.at(j).pa->coor)-dist3d(graph.at(i).pb->coor,graph.at(j).pb->coor));        
          if(dist < ps_dDist || fabs(dist-ps_dDist)<0.001){ //If passes distance threshold
            // cout<<endl<<graph.at(i).a->coor[0]<<" "<<graph.at(i).a->coor[1]<<" "<<graph.at(i).a->coor[2]<<" "<<graph.at(j).a->coor[0]<<" "<<graph.at(j).a->coor[1]<<" "<<graph.at(j).a->coor[2];
            // cout<<endl<<graph.at(i).b->coor[0]<<" "<<graph.at(i).b->coor[1]<<" "<<graph.at(i).b->coor[2]<<" "<<graph.at(j).b->coor[0]<<" "<<graph.at(j).b->coor[1]<<" "<<graph.at(j).b->coor[2];
            // cout<<endl<<"Dist: "<<dist<<endl<<i<<" "<<j<< " ConnID: "<< ConnID(i,j)<<" Connection!"<<endl;
            conn[ConnID(i,j)]=1;
            numEdges++; //create the new edge
          }
        }
      }else{
        if(graph.at(i).a->id != graph.at(j).a->id && graph.at(i).b->id != graph.at(j).b->id){
          dist=fabs(dist3d(graph.at(i).a->coor,graph.at(j).a->coor)-dist3d(graph.at(i).b->coor,graph.at(j).b->coor));        
          if(dist < dDist || fabs(dist-dDist)<0.001){ //If passes distance threshold
            // cout<<endl<<graph.at(i).a->coor[0]<<" "<<graph.at(i).a->coor[1]<<" "<<graph.at(i).a->coor[2]<<" "<<graph.at(j).a->coor[0]<<" "<<graph.at(j).a->coor[1]<<" "<<graph.at(j).a->coor[2];
            // cout<<endl<<graph.at(i).b->coor[0]<<" "<<graph.at(i).b->coor[1]<<" "<<graph.at(i).b->coor[2]<<" "<<graph.at(j).b->coor[0]<<" "<<graph.at(j).b->coor[1]<<" "<<graph.at(j).b->coor[2];
            // cout<<endl<<"Dist: "<<dist<<endl<<i<<" "<<j<< " ConnID: "<< ConnID(i,j)<<" Connection!"<<endl;
            conn[ConnID(i,j)]=1;
            numEdges++; //create the new edge
          }
        }
      }
    }
  }
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void bk(int cg, vector<node> &graph, bool* &conn){
  int c;
  int* all = new int[numNodes];
  compsub = new int[numNodes];

  for(c=0; c<numNodes; c++){
    all[c]=c;
  }

  // sortArray(all, numNodes, conn);

  Extend(all,0,numNodes,cg,graph,conn,0);

  delete[] compsub;
  compsub=NULL;
  
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void sortArray(int * &in, int nn, bool* &conn){
  vector<nodeI> sortn;
  int i=0;
  int c=0;

  //Sort array
  for(c=0;c<nn;c++){
    nodeI newn;
    newn.id=in[c];
    newn.neibrs=0;
    for(i=0;i<nn;i++){
      if(!conn[ConnID(in[c],i)]) newn.neibrs++;
    }
    sortn.push_back(newn);
    // cout<<endl<< newn.id<<" "<< newn.neibrs;
  }

  sort (sortn.begin(), sortn.end(), myfunction);

  vector<nodeI>::iterator pnit;
  c=0;
  for(pnit=sortn.begin(); pnit!=sortn.end(); ++pnit){
      in[c]=(*pnit).id;
      c++;
  }
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
bool myfunction (nodeI i,nodeI j) { return (i.neibrs<j.neibrs); }
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
 void Extend(int* old,int ne,int ce, int cg, vector<node> &graph, bool* &conn, int lev){
  int fixp;
  int newne,newce,j,count,pos,p,s,sel,loc,l;
  int* neww = new int[ce];
  int minnod=ce;
  int i=0;
  int nod=0;

  lev++;

  if(lev==1){ stopBk=0; }

  //Determine each count value and look for the one with least disconnections
  // cout<<endl<<"Entering Extend with "<<ce<<" candidates and "<<ne<<" nots lev "<<lev<<endl;
  while(i < ce && minnod != 0)
  {
    p = old[i]; // Id of current cand
    count = 0; // Reset its disconnections counter

    j = ne; // To count disconnections with candidates only
    // Count disconnections with candidates only
    while(j < ce)  //&& count < minnod
    {
      if(!conn[ConnID(p,old[j])]){
        count++; // Increment disconnection counter
        pos=j; // Set pivot as last candidate of the list
      }
      j++;
    }
    // Remark : first "old" set is All and n_not = 0, so fixed_point and s will first be the node with the most connections.
    if(count < minnod){ //If this node has less disconnections with previous one, set it as pivot
      fixp=p;
      minnod=count;
      if(i < ne){ // If its from not, set the pivot as the last candidate
        s=pos;
      }else{ // Else if its a candidate take it as pivot
        s=i;
        nod=1;
      }
    }
    i++;
  }

  // cout<<"Expanding "<<minnod+nod<<" times. bk "<<stopBk<<endl;
  // Tries to expand the current clique with each candidate starting with the one with least disconnections
  for(nod=minnod+nod;nod>0;nod--){

    // cout<<"ce "<<ce<<" ne "<<ne<<" lev "<<lev<<endl;
    // cout<<"New PIVOT "<<old[s]<<endl;
    // cout<<"stopbk "<<stopBk<<endl;

    if(stopBk==1) break;

    // for(int z=0; z<ce; z++){
      // if(z==ce) cout<<"| ";
      // cout<<old[z]<<" ";
    // }
    // cout<<endl;

    // Switch position of pivot with the first candidate
    p = old[s];
    old[s] = old[ne];
    old[ne] = p;
    sel = p;

    // for(int z=0; z<ce; z++){
      // if(z==ce) cout<<"| ";
      // cout<<old[z]<<" ";
    // }
    // cout<<endl;

    // Build the new "not" set based on connections between node sel and old "not" set
    newne = 0;
    i = 0;
    while(i<ne){
      // cout<<i<<" "<<old[i];
      if(conn[ConnID(sel,old[i])]){
        neww[newne++]=old[i];
        // cout<<" conn not";
      }
      // cout<<endl;
      i++;
    }

    // Build the new "candidates" set based on connections between node sel and old "candidates" set
    newce=newne;
    i=ne+1;
    while(i < ce){
      // cout<<i<<" "<<old[i];
      if(conn[ConnID(sel,old[i])]){
        neww[newce++]=old[i];
        // cout<<" conn can";
      }
      // cout<<endl;
      i++;
    }

    // Add node sel in clique
    compsub[c++]=sel;
    // cout<<"COMPSUB: ";
    // for(int g=0; g<c; g++){
      // cout<<compsub[g]<<" ";
    // }
    // cout<<endl;

      // cout<<"newce "<<newce<<" newne "<<newne<<" c "<<c<<" lev "<<lev<<endl;

    // cout<<"| ";
    // for(int g=0; g<newce; g++){
      // if(g==newne) cout<<"| ";
      // cout<<"["<<g<<"] "<<neww[g]<<" ";
    // }
    // cout<<endl;

    // if((newce+c)>=topN){ //Check if there is enough candidate in this extension to find a new TOP clique
      if(newce == 0){ // Print clique if both "not" and "candidates" sets are empty
        if(c >= Clique_threshold){
          AddNewClique(c,compsub,cg,graph);
          stopBk = 1;
          // cout<< nCliques<<endl;
          if(bkAll == 1 && (nCliques<maxCliques)){
            stopBk = 0;
          }
        }
      }else if(newne < newce){
        // Continue to expand clique if there are remaining candidates or "not" nodes
        Extend(neww,newne,newce,cg,graph,conn,lev);
      }      
    // }else{
      // int tmpint=newce+c;
      // cout<<"exiting max size "<<tmpint<<endl;
    // }


    // Collapse clique (Recursion will be complete and clique will be printed at this point)
    c--;

    // cout<<endl<<"COMPSUB: ";
    // for(int g=0; g<c; g++){
      // cout<<compsub[g]<<" ";
    // }

    // Place node sel in "not" set (Is the last recorded point of the printed clique)
    ne++;

    // cout<<endl<<"neOUT | ";
    // for(int g=0; g<ce; g++){
    //   if(g==ne) cout<<"| ceOUT | ";
    //   cout<<old[g]<<" ";
    // }

    // cout<<endl<<"Fixp "<<fixp<<" ne "<<ne<<" nod "<<nod<<endl;

    // Find new starting point (new node position s) for the recursion with no connection with node fixed_point (which was part of last clique and had the most connections)
    // May be able to find a new part of the clique through connections with candidates and disconnection with fixed_point
    if(nod > 1){
      s=ne;
      while(conn[ConnID(fixp,old[s])] && s < numNodes) s++;
    }

    if(lev>2) break;
  }
  lev--;
  // cout<<"exiting extend"<<endl;

  delete[] neww;
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void AddNewClique(int n, int* list, int cg, vector<node> &graph){
  vector<nodes>::iterator it;
  vector<float> la;
  vector<float> lb;

  nCliques++;
  // cout<<nCliques<<endl;

  Clique newClique;
  newClique.cg=cg;
  newClique.nbNodes=n;
  newClique.nbNodesM=0;
  newClique.tani=0.0;
  newClique.taniM=0.0;
  newClique.taniNorm=0.0;
  newClique.normNodes=0.0;
  newClique.normNodesRMSD=0.0;
  newClique.rmsd=0.0;
  newClique.cosdBu=0.0;
  newClique.nrg=0.0;

  for(int i=0; i<n; i++){ newClique.nodes.push_back(graph.at(list[i])); }
  cliques.push_back(newClique);

  //Get list of coords to get rotation matrix
  if(cg==-1){
    for(it=cliques.back().nodes.begin(); it!=cliques.back().nodes.end(); ++it){
      for(int i=0; i<3; i++){
        la.push_back((*it).ca->coor[i]);
        lb.push_back((*it).cb->coor[i]);
      }
    }
  }else if(cg==-3){
    for(it=cliques.back().nodes.begin(); it!=cliques.back().nodes.end(); ++it){
      for(int i=0; i<3; i++){
        la.push_back((*it).pa->coor[i]);
        lb.push_back((*it).pb->coor[i]);
      }
    }
  }else{
    for(it=cliques.back().nodes.begin(); it!=cliques.back().nodes.end(); ++it){
      for(int i=0; i<nb_of_probes; ++i){
        if((*it).a->pb[i]==1 && (*it).b->pb[i]==1){
          for(int j=0; j<3; j++){
            la.push_back((*it).a->coor[j]);
            lb.push_back((*it).b->coor[j]);
          }
        }
      }
    }
  }

  //Calculate rotation matrix
  cliques.back().mat_r=gsl_matrix_alloc(3,3);
  gsl_matrix_set_zero(cliques.back().mat_r);
  for(int i=0; i<3; i++){
    cliques.back().cen_a[i]=0.0;
    cliques.back().cen_b[i]=0.0;
  }
  cliques.back().det=calcRot(la,lb,cliques.back().cen_a,cliques.back().cen_b,cliques.back().mat_r);

  //Rotate mif 1 onto mif 2
  for(int v=0; v<mif1.size(); v++){
    for(int i=0; i<3; i++){
      mif1[v].ncoor[i]=cliques.back().cen_b[i];
      for(int j=0; j<3; j++){
        mif1[v].ncoor[i]+=(mif1[v].coor[j]-cliques.back().cen_a[j])*gsl_matrix_get(cliques.back().mat_r,i,j);
      }
    }
  }

  // Rotate ligand and calculate RMSD
  float ligRMSD=0.0;
  int ligRMSDc=0;
  if(rnc1.compare("")!=0 && rnc2.compare("")!=0 && lig1.size()>0 && lig2.size()>0){
    for(int v=0; v<lig1.size(); v++){
      float dist=0.0;
      for(int i=0; i<3; i++){
        lig1[v].ncoor[i]=cliques.back().cen_b[i];
        for(int j=0; j<3; j++){ lig1[v].ncoor[i]+=(lig1[v].coor[j]-cliques.back().cen_a[j])*gsl_matrix_get(cliques.back().mat_r,i,j); }
      }
      for(int w=0; w<lig2.size(); w++){
        if(lig2[w].atomn.compare(lig1[v].atomn)==0){
          ligRMSD+=pow(dist3d(lig1[v].ncoor,lig2[w].coor),2.0);
          // cout<<lig1[v].atomn<<" "<<lig2[w].atomn<<" "<<pow(dist3d(lig1[v].ncoor,lig2[w].coor),2.0)<<endl;
          ligRMSDc++;
          break;
        }
      }
    }
    if(ligRMSDc>0 && ligRMSDc==lig1.size() && lig1.size()==lig2.size()){
      ligRMSD=sqrt(ligRMSD/(float)ligRMSDc);
    }else{
      // cout<<"ligRMSDc "<<ligRMSDc<<" ligsize1 "<<lig1.size()<<" ligsize2 "<<lig2.size()<<endl;
      ligRMSD=0.0;
    }
  }
  cliques.back().ligRMSD=ligRMSD;
  // cout<<"LigRMSD "<<ligRMSD<<endl;

  //Calculate RMSD,vecCosineDiff
  float rmsd=0.0;
  float numBu=0.0;
  float d1Bu=0.0;
  float d2Bu=0.0;
  for(it=cliques.back().nodes.begin(); it!=cliques.back().nodes.end(); ++it){
    float ncoor[3];
    for(int i=0; i<3; i++){
      ncoor[i]=cliques.back().cen_b[i];
      for(int j=0; j<3; j++){
        if(cg==-1){
          ncoor[i]+=((*it).ca->coor[j]-cliques.back().cen_a[j])*gsl_matrix_get(cliques.back().mat_r,i,j);
        }else if(cg==-3){
          ncoor[i]+=((*it).pa->coor[j]-cliques.back().cen_a[j])*gsl_matrix_get(cliques.back().mat_r,i,j);
        }else{
          ncoor[i]+=((*it).a->coor[j]-cliques.back().cen_a[j])*gsl_matrix_get(cliques.back().mat_r,i,j);
        }
      }
    }
    if(cg==-1){
      rmsd+=pow(dist3d(ncoor,(*it).cb->coor),2.0);
    }else if(cg==-3){
      rmsd+=pow(dist3d(ncoor,(*it).pb->coor),2.0);
    }else{
      rmsd+=pow(dist3d(ncoor,(*it).b->coor),2.0);

      // cout<<(*it).a->bu<<" "<<(*it).b->bu<<endl;
      numBu+=(*it).a->bu*(*it).b->bu;
      d1Bu+=pow((*it).a->bu,2);
      d2Bu+=pow((*it).b->bu,2);

      for(int i=0; i<nb_of_probes; ++i){
        if((*it).a->pb[i]==1 && (*it).b->pb[i]==1) cliques.back().nrg+=(*it).a->nrg[i]+(*it).b->nrg[i];
      }

      cliques.back().nbNodesM+=(*it).nbi;

      // printf("A %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",(*it).a->nrg[0],(*it).a->ang[0],(*it).a->nrg[1],(*it).a->ang[1],(*it).a->nrg[2],(*it).a->ang[2],(*it).a->nrg[3],(*it).a->ang[3],(*it).a->nrg[4],(*it).a->ang[4],(*it).a->nrg[5],(*it).a->ang[5]);
      // printf("B %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",(*it).b->nrg[0],(*it).b->ang[0],(*it).b->nrg[1],(*it).b->ang[1],(*it).b->nrg[2],(*it).b->ang[2],(*it).b->nrg[3],(*it).b->ang[3],(*it).b->nrg[4],(*it).b->ang[4],(*it).b->nrg[5],(*it).b->ang[5]);
      // printf("A %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",(*it).a->nrg[0],(*it).a->nrg[1],(*it).a->nrg[2],(*it).a->nrg[3],(*it).a->nrg[4],(*it).a->nrg[5]);
      // printf("B %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",(*it).b->nrg[0],(*it).b->nrg[1],(*it).b->nrg[2],(*it).b->nrg[3],(*it).b->nrg[4],(*it).b->nrg[5]);
      // cout<<"cosd "<<(*it).cosd<<endl;
      cliques.back().normNodes+=(*it).cosd;
    }
  }
  cliques.back().cosdBu=numBu/(sqrt(d1Bu)*sqrt(d2Bu));
  cliques.back().rmsd=sqrt(rmsd/(float)cliques.back().nbNodes);
  cliques.back().normNodesRMSD=cliques.back().normNodes/cliques.back().rmsd;
  
  // cout<<"Finding corresponding vertexes..."<<endl;
  // float dist=0.0;
  // for(int u=0; u<mif1.size(); u++){
  //   if(mif1[u].grid[cg2]!=1) continue;
  //   for(it=cliques.back().nodes.begin(); it!=cliques.back().nodes.end(); ++it){
  //     dist=dist3d(mif1[u].ncoor,(*it).a.coor);
  //   }
  //   for(int v=0; v<mif2.size(); v++){
  //     if(mif2[v].grid[cg2]!=1) continue;
  //     dist=dist3d(mif1[u].ncoor,mif2[v].coor);  
  //     if(dist < dDist || fabs(dist-dDist)<0.001){ //If passes distance threshold
  //       for(int i=0; i<nb_of_probes; i++){
  //         if(mif1[u].pb[i]==1 && mif2[v].pb[i]==1){
  //           cout<<i<<" - "<<mif1[u].ncoor[0]<<" "<<mif1[u].ncoor[1]<<" "<<mif1[u].ncoor[2]<<" "<<mif2[v].coor[0]<<" "<<mif2[v].coor[1]<<" "<<mif2[v].coor[2]<<" - "<<mif1[u].pb[i]<<" "<<mif2[v].pb[i]<<endl;  
  //           mif1[u].m[i]=1;
  //           mif2[v].m[i]=1;
  //         }
  //       }
  //     }
  //   }
  // }

  if(cg==-1){
    cliques.back().tani=(float)n/((float)caSize1+(float)caSize2-(float)n);
  }else if(cg==-3){
    cliques.back().tani=(float)n/((float)cliques.back().rmsd);
  }else{
    cliques.back().tani=(float)n/((float)ss1[cg]+(float)ss2[cg]-(float)n);
    cliques.back().taniM=(float)cliques.back().nbNodesM/((float)ss1m[cg]+(float)ss2m[cg]-(float)cliques.back().nbNodesM);
    cliques.back().taniNorm=cliques.back().normNodes/((float)ss1[cg]+(float)ss2[cg]-cliques.back().normNodes);
  }

  // cout<<"ROTATION MAT:";
  // for(int i=0; i<3; i++){
    // for(int j=0; j<3; j++){
      // cout<<" "<<gsl_matrix_get(cliques.back().mat_r,i,j);
    // }
  // }
  // cout<<endl<<"cen_a: "<<cliques.back().cen_a[0]<<" "<<cliques.back().cen_a[1]<<" "<<cliques.back().cen_a[2];
  // cout<<endl<<"cen_b: "<<cliques.back().cen_b[0]<<" "<<cliques.back().cen_b[1]<<" "<<cliques.back().cen_b[2]<<endl;
  
  if(cliques.back().taniNorm>topT){
    cout<<"NEW TOP CLIQUE CG "<<cg<<" cosdBu "<<cliques.back().cosdBu<<" nbNodes "<<cliques.back().nbNodes<<" cInts "<<cliques.back().nbNodesM<<" nrg "<<cliques.back().nrg<<" normNodes "<<cliques.back().normNodes<<" normNodesRMSD "<<cliques.back().normNodesRMSD<<" tani "<<cliques.back().tani<<" taniM "<<cliques.back().taniM<<" taniNormNodes "<<cliques.back().taniNorm<<" RMSD "<<cliques.back().rmsd<<" ligRMSD "<<cliques.back().ligRMSD<<endl;
    topT=cliques.back().taniNorm;
    topN=cliques.back().nbNodes;
    topCliques[cg]=cliques.size()-1;
  }else{
    cout<<"CLIQUE CG "<<cg<<" cosdBu "<<cliques.back().cosdBu<<" nbNodes "<<cliques.back().nbNodes<<" cInts "<<cliques.back().nbNodesM<<" nrg "<<cliques.back().nrg<<" normNodes "<<cliques.back().normNodes<<" normNodesRMSD "<<cliques.back().normNodesRMSD<<" tani "<<cliques.back().tani<<" taniM "<<cliques.back().taniM<<" taniNormNodes "<<cliques.back().taniNorm<<" RMSD "<<cliques.back().rmsd<<" ligRMSD "<<cliques.back().ligRMSD<<endl;
  }

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void printNodes(){
  
  vector<nodes>::iterator it;
  char suffix[50];

  if(pairwiseF.compare("")==0){
    if(wrfn==1){ //Add similarity score to filename
      if(rnc1.compare("")!=0 && rnc2.compare("")!=0){
        sprintf(suffix,"_%d_%5.4f_%5.4f",cliques[topCliques[steps.back()]].nbNodes,cliques[topCliques[steps.back()]].tani,cliques[topCliques[steps.back()]].ligRMSD);
      }else{
        sprintf(suffix,"_%d_%5.4f",cliques[topCliques[steps.back()]].nbNodes,cliques[topCliques[steps.back()]].tani);
      }
      strcat(out_file,suffix);
    }
    strcat(out_file,".isomif");
    open_file_ptr(&fpout,out_file,1);
  }


  fprintf(fpout,"%s",outH);
  for(int cs=0; cs<steps.size(); cs++){
    int istart=0;
    int iend=cliques.size();
    if(wc==0){
      istart=topCliques[steps[cs]];
      iend=topCliques[steps[cs]]+1;
    }
    for(int i=istart; i<iend; i++){
      if(steps[cs]==-2){
        fprintf(fpout,"REMARK CLIQUE CG %d NODES %d TANI %5.3f SS1 %d SS2 %d\n",cliques[i].cg,cliques[i].nbNodes,cliques[i].tani,ss1[cg2],ss2[cg2]);
        if(emptOut!=1){
          for(int j=0; j<cliques[i].va.size(); j++){
            fprintf(fpout, "A %8.3f %8.3f %8.3f %d %d %d %d %d %d\n",cliques[i].va[j].ncoor[0],cliques[i].va[j].ncoor[1],cliques[i].va[j].ncoor[2],cliques[i].va[j].m[0],cliques[i].va[j].m[1],cliques[i].va[j].m[2],cliques[i].va[j].m[3],cliques[i].va[j].m[4],cliques[i].va[j].m[5]);
          }
          for(int j=0; j<cliques[i].vb.size(); j++){
            fprintf(fpout, "B %8.3f %8.3f %8.3f %d %d %d %d %d %d\n",cliques[i].vb[j].coor[0],cliques[i].vb[j].coor[1],cliques[i].vb[j].coor[2],cliques[i].vb[j].m[0],cliques[i].vb[j].m[1],cliques[i].vb[j].m[2],cliques[i].vb[j].m[3],cliques[i].vb[j].m[4],cliques[i].vb[j].m[5]);
          }
        }
      }else if(steps[cs]==-1){
        fprintf(fpout,"REMARK CLIQUE CG %d NODES %d TANI %5.3f SS1 %d SS2 %d\n",cliques[i].cg,cliques[i].nbNodes,cliques[i].tani,(int)prot1.size(),(int)prot2.size());
        if(emptOut!=1){
          for(it=cliques[i].nodes.begin(); it!=cliques[i].nodes.end(); ++it){
            fprintf(fpout, "%3s %4d %4s %5d %s %8.3f %8.3f %8.3f %3s %4d %4s %5d %s %8.3f %8.3f %8.3f\n",(*it).ca->resn.c_str(),(*it).ca->resnb,(*it).ca->atomn.c_str(),(*it).ca->atomnb,(*it).ca->chain.c_str(),(*it).ca->coor[0],(*it).ca->coor[1],(*it).ca->coor[2],(*it).cb->resn.c_str(),(*it).cb->resnb,(*it).cb->atomn.c_str(),(*it).cb->atomnb,(*it).cb->chain.c_str(),(*it).cb->coor[0],(*it).cb->coor[1],(*it).cb->coor[2]);
          }
        }
      }else if(steps[cs]==-3){
        fprintf(fpout,"REMARK CLIQUE CG %d NODES %d TANI %5.3f SS1 %d SS2 %d\n",cliques[i].cg,cliques[i].nbNodes,cliques[i].tani,(int)pseudoL1.size(),(int)pseudoL2.size());
        if(emptOut!=1){
          for(it=cliques[i].nodes.begin(); it!=cliques[i].nodes.end(); ++it){
            fprintf(fpout, "%3s %8.3f %8.3f %8.3f %3s %8.3f %8.3f %8.3f\n",(*it).pa->type.c_str(),(*it).pa->coor[0],(*it).pa->coor[1],(*it).pa->coor[2],(*it).pb->type.c_str(),(*it).pb->coor[0],(*it).pb->coor[1],(*it).pb->coor[2]);
          }
        }
      }else{
        fprintf(fpout,"REMARK CLIQUE CG %d NODES %d NODESM %d NORMNODES %6.3f NORMNODESRMSD %6.3f TANI %5.3f TANIM %5.3f TANINORM %5.3f RMSD %5.3f NRG %.3f SS1 %d SS2 %d SS1M %d SS2M %d LIGRMSD %5.3f\n",cliques[i].cg,cliques[i].nbNodes,cliques[i].nbNodesM,cliques[i].normNodes,cliques[i].normNodesRMSD,cliques[i].tani,cliques[i].taniM,cliques[i].taniNorm,cliques[i].rmsd,cliques[i].nrg,ss1[cliques[i].cg],ss2[cliques[i].cg],ss1m[cliques[i].cg],ss2m[cliques[i].cg],cliques[i].ligRMSD);
        if(emptOut!=1){
          for(it=cliques[i].nodes.begin(); it!=cliques[i].nodes.end(); ++it){
            for(int j=0; j<nb_of_probes; ++j){
              if((*it).a->pb[j]==1 && (*it).b->pb[j]==1){
                fprintf(fpout, "%d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",j,(*it).a->coor[0],(*it).a->coor[1],(*it).a->coor[2],(*it).b->coor[0],(*it).b->coor[1],(*it).b->coor[2]);
              }
            }
          }
        }
      }
      if(emptOut!=1){
        fprintf(fpout,"REMARK ROTMAT ");
        for(int m=0; m<3; m++) {
          for(int n=0; n<3; n++){
            fprintf(fpout," %9.4f",gsl_matrix_get(cliques[i].mat_r,m,n));
          }
        }
        fprintf(fpout,"\nREMARK CENTRES %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",cliques[i].cen_a[0],cliques[i].cen_a[1],cliques[i].cen_a[2],cliques[i].cen_b[0],cliques[i].cen_b[1],cliques[i].cen_b[2]);
        fprintf(fpout,"REMARK DET %g\n",cliques[i].det);
      }
    }
    fprintf(fpout,"REMARK END\n");
  }
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void clearStep(int cg){
  vector<nodes>::iterator it;
  vector<float> la;
  vector<float> lb;

  if(topT==-1.0){
    Clique newClique;
    newClique.cg=cg;
    newClique.nbNodes=0;
    newClique.normNodes=0.0;
    newClique.normNodesRMSD=0.0;
    newClique.taniNorm=0.0;
    newClique.tani=0.0;
    newClique.ligRMSD=0.0;
    newClique.ligRMSD=0.0;
    newClique.mat_r=gsl_matrix_alloc(3,3);
    for(int i=0; i<3; i++) {
      newClique.cen_a[i]=0.0;
      newClique.cen_b[i]=0.0;
    }
    gsl_matrix_set_zero(newClique.mat_r);
    cliques.push_back(newClique);
    topCliques[cg]=cliques.size()-1;
  }

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
double calcRot(vector<float> lista, vector<float> listb, float* cen_a, float* cen_b, gsl_matrix *mat_r){
  int i,j,k;
  
  // for(j=0; j<lista.size(); j+=3){
  //   cout<<lista.at(j)<<" "<<lista.at(j+1)<<" "<<lista.at(j+2)<<" | "<<listb.at(j)<<" "<<listb.at(j+1)<<" "<<listb.at(j+2)<<endl;
  // }

  //Get average X, Y and Z
  for(i=0; i<3; i++){
    cen_a[i] = 0.0;
    cen_b[i] = 0.0;
    for(j=0; j<lista.size(); j+=3){
      cen_a[i] += lista.at(j+i);
      cen_b[i] += listb.at(j+i);
    }
    cen_a[i] /= lista.size()/3;
    cen_b[i] /= lista.size()/3;

    // cout<<i<<" cen_a "<<cen_a[i]<<" cen_b "<<cen_b[i]<<endl;
  }

  // Create Co-Variance Matrix
  float val=0.0;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      val = 0.0;
      for(k=0; k<lista.size(); k+=3){
        val += (lista.at(k+i) - cen_a[i])*(listb.at(k+j) - cen_b[j]);
      }
      gsl_matrix_set(mat_r,i,j,val);
      // printf("%10.4f",gsl_matrix_get(mat_r,i,j));
    }
    // printf("\n");
  }
  return(SupSVD(mat_r));
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
double SupSVD(gsl_matrix *mat_u){
  int i,j,k,m;
  double val;
  double eps = 0.0001;

  gsl_matrix *mat_ut = gsl_matrix_alloc(3,3); 
  gsl_matrix *mat_v  = gsl_matrix_alloc(3,3);
  gsl_vector *vec_s  = gsl_vector_alloc(3);
  gsl_vector *vec_w  = gsl_vector_alloc(3);

  
  gsl_linalg_SV_decomp(mat_u,mat_v,vec_s,vec_w);
  /*
    This function factorizes the M-by-N matrix A (mat_u) into the singular 
    value  decomposition A = USV^T. On output the matrix A is replaced by U. 
    The diagonal elements of the singular value matrix S are stored in the 
    vector S (vec_s). The singular values are non-negative and form a 
    nonincreasing sequence from S_1 to S_N. The matrix V (mat_v) contains the 
    elements of V in untransposed form. To form the product USV^T it is 
    necessary to take the transpose of V. A workspace of length N is required 
    in work (vec_w). This routine uses the Golub-Reinsch SVD algorithm. 
  */


  //Store the transpose of mat_u in mat_ut
  gsl_matrix_transpose_memcpy(mat_ut,mat_u);
  
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++){
      val=0.0;
      for(k=0;k<3;k++){
        val += gsl_matrix_get(mat_v,i,k)*gsl_matrix_get(mat_ut,k,j);
      }
      if(val*val < eps) val = 0.0;
      gsl_matrix_set(mat_u,i,j,val);      
    }
  }

  // cout<<endl<<"mat_u"<<endl;
  // for(int i=0; i<3; i++) {
  //   for(int j=0; j<3; j++){
  //     printf("%7.3f ",gsl_matrix_get(mat_u,i,j));
  //   }
  //   printf("\n");
  // }

  /* Add by DD 30-08-2013
  * 
  * Fixed problem with mirror effect. If determinant of V*U^t is -1
  * 
  */

  double DetU;
  DetU = gsl_matrix_Det3D(mat_u);
  if (DetU < 0.0) { 
    gsl_matrix_set(mat_v,0,2,-1.0 * gsl_matrix_get(mat_v,0,2));
    gsl_matrix_set(mat_v,1,2,-1.0 * gsl_matrix_get(mat_v,1,2));
    gsl_matrix_set(mat_v,2,2,-1.0 * gsl_matrix_get(mat_v,2,2));
    for(i=0;i<3;i++) {
      for(j=0;j<3;j++){
        m=3*i+j;
        val=0.0;
        for(k=0;k<3;k++){
          val += gsl_matrix_get(mat_v,i,k)*gsl_matrix_get(mat_ut,k,j);
        }
        if(val*val < eps) val = 0.0;
        gsl_matrix_set(mat_u,i,j,val);      
      }
    }
  }

  gsl_matrix_free(mat_ut);
  gsl_matrix_free(mat_v);
  gsl_vector_free(vec_s);
  gsl_vector_free(vec_w);
  
  return gsl_matrix_Det3D(mat_u);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
double gsl_matrix_Det3D(gsl_matrix *M){
  int i,j;
  double det;

  //  guide to indexes: 0=x, 1=y, 2=z
  //  00 01 02
  //  10 11 12
  //  20 21 22

  det = gsl_matrix_get(M,0,0)*gsl_matrix_get(M,1,1)*gsl_matrix_get(M,2,2)
    + gsl_matrix_get(M,0,1)*gsl_matrix_get(M,1,2)*gsl_matrix_get(M,2,0)
    + gsl_matrix_get(M,0,2)*gsl_matrix_get(M,2,1)*gsl_matrix_get(M,1,0)
    - gsl_matrix_get(M,0,2)*gsl_matrix_get(M,1,1)*gsl_matrix_get(M,2,0)
    - gsl_matrix_get(M,1,2)*gsl_matrix_get(M,2,1)*gsl_matrix_get(M,0,0)
    - gsl_matrix_get(M,2,2)*gsl_matrix_get(M,1,0)*gsl_matrix_get(M,0,1);

  return det;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void createNodes(int cg, vector<node> &graph, int s){

  float dist;
  if(cg==-1){
    for(int i=0; i<prot1.size(); i++){
      if(prot1.at(i).atomn.compare("CA")!=0 || prot1.at(i).bs!=1) continue; //needs to be a CA and in the binding site
      for(int j=0; j<prot2.size(); j++){
        if(prot2.at(j).atomn.compare("CA")!=0 || prot2.at(j).bs!=1) continue; //needs to be a CA and in the binding site
        // cout<< prot1.at(i).resn<<" "<< prot1.at(i).atomn<< " with "<< prot2.at(j).resn<<" "<< prot2.at(j).atomn;
        if(jtt[nam2n(prot1.at(i).resn)][nam2n(prot2.at(j).resn)]==1){ //They are similar based on the JTT matrix and JTT threshold
          node newNode;
          newNode.ca=&prot1.at(i);
          newNode.cb=&prot2.at(j);
          graph.push_back(newNode);
        }else{
          // cout<< " 0"<<endl;
        }
      }
    }
  }else if(cg==-3){
    for(int i=0; i<pseudoL1.size(); i++){
      for(int j=0; j<pseudoL2.size(); j++){
        if(samePseudo(pseudoL1.at(i),pseudoL2.at(j))==1){ //They are similar based on the JTT matrix and JTT threshold
          node newNode;
          newNode.pa=&pseudoL1.at(i);
          newNode.pb=&pseudoL2.at(j);
          graph.push_back(newNode);
        }
      }
    }
  }else{
    for(int i=0; i<mif1.size(); i++){
      if(mif1.at(i).grid[cg]!=1) continue;

      for(int j=0; j<mif2.size(); j++){
        if(mif2.at(j).grid[cg]!=1) continue;

          int flag=0;
          int cints=0;
          for(int pb=0; pb<nb_of_probes; pb++){
            if(mif1.at(i).pb[pb]==1 && mif2.at(j).pb[pb]==1) flag++;
            if(mif1.at(i).pb[pb]==mif2.at(j).pb[pb]) cints++;
          }

          float num=0.0;
          float d1=0.0;
          float d2=0.0;
          for(int pb=0; pb<nb_of_probes; ++pb){
            num+=(mif1.at(i).nrg[pb]*mif2.at(j).nrg[pb]);
            d1+=pow(mif1.at(i).nrg[pb],2);
            d2+=pow(mif2.at(j).nrg[pb],2);
          }
          float cosd=num/(sqrt(d1)*sqrt(d2));
          // cout<<"cosd: "<<cosd<<endl;
          // if(cosd<cosdT) continue;

          // cout<<cints<<" "<<commonInt<<" flag "<<flag<<endl;
          if(cints>=commonInt && flag!=0){

            // printf("A %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",mif1.at(i).nrg[0],mif1.at(i).ang[0],mif1.at(i).nrg[1],mif1.at(i).ang[1],mif1.at(i).nrg[2],mif1.at(i).ang[2],mif1.at(i).nrg[3],mif1.at(i).ang[3],mif1.at(i).nrg[4],mif1.at(i).ang[4],mif1.at(i).nrg[5],mif1.at(i).ang[5]);
            // printf("B %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",mif2.at(j).nrg[0],mif2.at(j).ang[0],mif2.at(j).nrg[1],mif2.at(j).ang[1],mif2.at(j).nrg[2],mif2.at(j).ang[2],mif2.at(j).nrg[3],mif2.at(j).ang[3],mif2.at(j).nrg[4],mif2.at(j).ang[4],mif2.at(j).nrg[5],mif2.at(j).ang[5]);
            // cout<<"cosd "<<cosd<<endl<<endl;

            if(s==1){
              // cout<<mif1.at(i).ncoor[0]<<" "<<mif1.at(i).ncoor[1]<<" "<<mif1.at(i).ncoor[2]<<" | "<<mif2.at(j).coor[0]<<" "<<mif2.at(j).coor[1]<<" "<<mif2.at(j).coor[2];
            
              //Verify if superimposed vrtx is within the distance threshold of vrtx of mif 2
              dist=dist3d(mif1.at(i).ncoor,mif2.at(j).coor);
              if(dist>neibr_dDist){
                // cout<<" dist "<<dist<<" >thresh"<<endl;
                continue;
              }else{
                // cout<<" dist "<<dist<<" OK!"<<endl;
              }    
            }
            // if(cosd<0.50) continue;
            // cout<<mif1.at(i).ncoor[0]<<" "<<mif1.at(i).ncoor[1]<<" "<<mif1.at(i).ncoor[2]<<" | "<<mif1.at(i).pb[0]<<" "<<mif1.at(i).pb[1]<<" "<<mif1.at(i).pb[2]<<" "<<mif1.at(i).pb[3]<<" "<<mif1.at(i).pb[4]<<" "<<mif1.at(i).pb[5]<<" <- SAME -> ";
            // cout<<mif2.at(j).coor[0]<<" "<<mif2.at(j).coor[1]<<" "<<mif2.at(j).coor[2]<<" | "<<mif2.at(j).pb[0]<<" "<<mif2.at(j).pb[1]<<" "<<mif2.at(j).pb[2]<<" "<<mif2.at(j).pb[3]<<" "<<mif2.at(j).pb[4]<<" "<<mif2.at(j).pb[5]<<endl;
            // printf("A %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",mif1.at(i).nrg[0],mif1.at(i).ang[0],mif1.at(i).nrg[1],mif1.at(i).ang[1],mif1.at(i).nrg[2],mif1.at(i).ang[2],mif1.at(i).nrg[3],mif1.at(i).ang[3],mif1.at(i).nrg[4],mif1.at(i).ang[4],mif1.at(i).nrg[5],mif1.at(i).ang[5]);
            // printf("B %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",mif2.at(j).nrg[0],mif2.at(j).ang[0],mif2.at(j).nrg[1],mif2.at(j).ang[1],mif2.at(j).nrg[2],mif2.at(j).ang[2],mif2.at(j).nrg[3],mif2.at(j).ang[3],mif2.at(j).nrg[4],mif2.at(j).ang[4],mif2.at(j).nrg[5],mif2.at(j).ang[5]);
            // cout<<"cosd: "<<cosd<<endl;

            node newNode;
            newNode.a=&mif1.at(i);
            newNode.b=&mif2.at(j);
            newNode.cosd=cosd;
            newNode.nbi=flag;
            graph.push_back(newNode);
          }
      }
    }
  }

  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int samePseudo(pseudoC& ps1,pseudoC& ps2){
  if(ps1.type.compare("don")==0 && ps2.type.compare("don")==0){
    return(1);
  }else if(ps1.type.compare("acc")==0 && ps2.type.compare("acc")==0){
    return(1);
  }else if(ps1.type.compare("doa")==0 && ps2.type.compare("acc")==0){
    return(1);
  }else if(ps1.type.compare("doa")==0 && ps2.type.compare("don")==0){
    return(1);
  }else if(ps1.type.compare("acc")==0 && ps2.type.compare("doa")==0){
    return(1);
  }else if(ps1.type.compare("don")==0 && ps2.type.compare("doa")==0){
    return(1);
  }else if(ps1.type.compare("doa")==0 && ps2.type.compare("doa")==0){
    return(1);
  }else if(ps1.type.compare("hyd")==0 && ps2.type.compare("hyd")==0){
    return(1);
  }else if(ps1.type.compare("arm")==0 && ps2.type.compare("arm")==0){
    return(1);
  }
  return(0);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void getPairwise(){
  string empt ("");
  string line;
  ifstream infile(pairwiseF.c_str());
  while(getline(infile,line)){
    if(line.compare("")!=0){
      istringstream ss(line);
      istream_iterator<string> begin(ss), end;
      vector<string> vec(begin, end);
      string mif1F="";
      string mif2F="";
      string rnc1id="";
      string rnc2id="";
      for(int i=0; i<vec.size(); i+=2){
        if(vec[i].compare("-p1")==0) mif1F=vec[i+1];
        if(vec[i].compare("-p2")==0) mif2F=vec[i+1];
        if(vec[i].compare("-l1")==0) rnc1id=vec[i+1];
        if(vec[i].compare("-l2")==0) rnc2id=vec[i+1];
      }
      pwRun npw;
      npw.mif1=mif1F;
      npw.mif2=mif2F;
      npw.rnc1=rnc1id;
      npw.rnc2=rnc2id;
      pw.push_back(npw);
    }
  }
  //Store pairwise filename
  const char* prefix3; //Pairwise prefix
  string pwpre;
  int i;
  for(i=pairwiseF.length()-1; i>=0; i--){
      if(pairwiseF.at(i)=='/'){
       pwpre=pairwiseF.substr(i+1,pairwiseF.length()-i-1);
       break;
     }
   }
  if(pwpre.compare(empt)==0){
    pwpre=pairwiseF.substr(0,pairwiseF.length()-i-1);
  }
  prefix3 = pwpre.c_str();

  if(!strcmp(outbase,"")){ //If no outbase name is defined, create one
    sprintf(out_file,"./pw_%s",prefix3);
  }else{
    sprintf(out_file,"%spw_%s",outbase,prefix3);
  }
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int createVrtxVec(string mifFile, vector<vertex>& p, vector<atom>& a, int* ss, int* ssm, int &caSize, vector<pseudoC>& pl, string rnc, vector<atom>& llist){
  string line;
  vertex* nvrtx;
  float x,y,z;
  string resn;
  int resnb,mif,bs,bu;
  int atmC=0;
  int pcC=0;
  string atomn;
  string pseudo;
  string thisresnumc;
  int atomnb;
  string chain;
  string dump;
  // int pb0,pb1,pb2,pb3,pb4,pb5,gr0,gr1,gr2,gr3;
  // float nrg0,nrg1,nrg2,nrg3,nrg4,nrg5,ang0,ang1,ang2,ang3,ang4,ang5;
  // int e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15,e16,e17,e18,e19;

  ifstream infile(mifFile.c_str());
  while(getline(infile,line)){
    if(line.compare("")==0){
      continue;
    }else if(line.compare(0,11,"#nbOfProbes")==0){
      stringstream test(line);
      test >> dump >> nb_of_probes;
    }else if(line.compare(0,5,"#ATOM")==0){
      stringstream test(line);
      test >> dump >> resn >> resnb >> atomn >> atomnb >> chain >> x >> y >> z >> mif >> bs;
      atom natom;
      natom.atomn=atomn;
      natom.atomnb=atomnb;
      natom.resn=resn;
      natom.resnb=resnb;
      natom.chain=chain;
      natom.coor[0]=x;
      natom.coor[1]=y;
      natom.coor[2]=z;
      natom.id=atmC;
      natom.mif=mif;
      natom.bs=bs;
      stringstream sss;
      sss << resnb;
      thisresnumc = resn + sss.str() + chain;
      if(rnc.compare(thisresnumc)==0){
        size_t found;
        found = atomn.find("H");
        if(found==string::npos) llist.push_back(natom);
      }
      if(atomn.compare("CA")==0 && bs==1) caSize++;
      a.push_back(natom);
      atmC++;
    }else if(line.compare(0,7,"#PSEUDO")==0){
      stringstream test(line);
      test >> dump >> pseudo >> x >> y >> z;
      pseudoC npseudo;
      npseudo.type=pseudo;
      npseudo.coor[0]=x;
      npseudo.coor[1]=y;
      npseudo.coor[2]=z;
      npseudo.id=pcC;
      pl.push_back(npseudo);
      pcC++;
    }else if(line.compare(0,1,"#")!=0){
      istringstream iss(line);
      vector<string> tokens;
      copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(tokens));

      // for(int i=1; i<tokens.size(); i++){
      //   eps[coln[i-1]+"_"+rown]=atoi(tokens[i].c_str());
      //   eps[rown+"_"+coln[i-1]]=atoi(tokens[i].c_str());
      // }

      // stringstream test(line);
      // test >> x >> y >> z >>pb0 >>nrg0 >>ang0 >>pb1 >>nrg1 >>ang1 >>pb2 >>nrg2 >>ang2 >>pb3 >>nrg3 >>ang3 >>pb4 >>nrg4 >>ang4 >>pb5 >>nrg5 >>ang5 >>gr0 >>gr1 >>gr2 >>gr3 >>bu;
      vertex nvrtx;
      nvrtx.coor[0]=atof(tokens[0].c_str());
      nvrtx.coor[1]=atof(tokens[1].c_str());
      nvrtx.coor[2]=atof(tokens[2].c_str());

      int pb=0;
      int intf=0;
      // cout<<line<<endl;
      // cout<<"token size "<<tokens.size()<<endl;
      while(pb<nb_of_probes){
        int pbi=(pb*3)+3;
        // cout<<pb<<" "<<pbi<<" "<<atoi(tokens[pbi].c_str())<<" "<<atof(tokens[pbi+1].c_str())<<" "<<atof(tokens[pbi+2].c_str())<<endl;
        if(atoi(tokens[pbi].c_str())==1) intf=1; //flag for search space
        for(int k=0; k<4; k++){
          if(nvrtx.grid[k]==1){
            if(atoi(tokens[pbi].c_str())==1) ssm[k]++; //increment probe search space
          }
        }
        nvrtx.pb.push_back(atoi(tokens[pbi].c_str()));
        nvrtx.nrg.push_back(atof(tokens[pbi+1].c_str()));
        nvrtx.ang.push_back(atof(tokens[pbi+2].c_str()));
        nvrtx.m.push_back(0);
        pb++;
      }

      nvrtx.grid[0]=atoi(tokens[tokens.size()-5].c_str());
      nvrtx.grid[1]=atoi(tokens[tokens.size()-4].c_str());
      nvrtx.grid[2]=atoi(tokens[tokens.size()-3].c_str());
      nvrtx.grid[3]=atoi(tokens[tokens.size()-2].c_str());
      nvrtx.cg.push_back(0);
      nvrtx.cg.push_back(0);
      nvrtx.cg.push_back(0);
      nvrtx.cg.push_back(0);
      nvrtx.id=totalVrtx; //Unique id for both clefts
      nvrtx.bu=atoi(tokens[tokens.size()-1].c_str());

      //Calculate search space
      if(intf==1){
        for(int k=0; k<4; k++){
          if(nvrtx.grid[k]==1){
            ss[k]++;
          }
        }
      }

      p.push_back(nvrtx);
      totalVrtx++;
    }
  }
  return(0);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
// bool compareSim(node a, node b) {
//   return a.similarity < b.similarity;
// }
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int open_file_ptr(FILE** file_ptr, char* filename, int open_type){
  if(open_type==1){
    (*file_ptr) = fopen(filename, "w");
  }else if(open_type==2){
    (*file_ptr) = fopen(filename, "a");
  }else if(open_type==3){
    (*file_ptr) = fopen(filename, "r");
  }
  if((*file_ptr) == NULL){
    printf("Can't open file %s\n",filename);
    return(24);
  }
  return(0);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
float dist3d(float a[], float b[]){
  float dist=0.0;
  int i;

  for(i=0;i<3;i++){
    dist += (a[i]-b[i])*(a[i]-b[i]);
  }
  dist = sqrt(dist);

  return dist;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int read_commandline(int argc, char *argv[]){

  char  usage[3000];
  char  tmp_line[250];
  int nb_arg;

  // assignment of default values to optional parameters
  rnc1="";
  rnc2="";

  strcpy(usage,"\n!---   IsoMIF   ---!\nWelcome.Bienvenue.\n");
  strcat(usage,"\nObligatory Arguments:\n");
  sprintf(tmp_line,"-p1          : \t Isomif energy filename of Protein 1\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-p2          : \t Isomif energy filename of Protein 2\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-pp          : \t File with highthroughput comparions\n");
  strcat(usage,tmp_line);
  strcat(usage,"\nOptional Arguments:\n");
  sprintf(tmp_line,"-h          : \t Print help menu\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-d          : \t Node Distance Threshold [default %4.2f A]\n",dDist);
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-dca        : \t C-alpha Distance Threshold [default %4.2f A]\n",ca_dDist);
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-dps        : \t Pseudocenter Distance Threshold [default %4.2f A]\n",ps_dDist);
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-dcn        : \t C-alpha-Node Distance Threshold [default %4.2f A]\n",neibr_dDist);
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-s          : \t Get all C-alpha cliques (1=yes, 0=no) [default %d]\n",bkAll);
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-c          : \t Coarse-grain step sequence -2=basic superimposition, -1=C-alpha, 0=(2.0 A), 1=(1.5 A), 2=(1.0 A) or 3=(0.5 A)\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-m          : \t Max nodes number per graph [default %d]\n",maxNodes);
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-o          : \t Output Directory Path [Default is generated using the Mif filenames]\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-l1         : \t Resnumc of ligand 1\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-l2         : \t Resnumc of ligand 2\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-q          : \t List of atoms from p1 and p2 to superimpose for step -2\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-w          : \t Print similarity in filename\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-e          : \t Empty result file\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-a          : \t Max Cliques\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-n          : \t Number of common interactions for a vector to be similar\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-j          : \t JTT matrix rank threhsold [default 5]\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-wc         : \t write all the cliques found and their info without the nodes [default 0]\n");
  strcat(usage,tmp_line);
  sprintf(tmp_line,"-pr         : \t print details on console\n");
  strcat(usage,tmp_line);


  //print command line in cmdLine
  strcpy(cmdLine,argv[0]);
  strcat(cmdLine," ");
  for(nb_arg=1; nb_arg<argc; nb_arg++){
    strcat(cmdLine,argv[nb_arg]);
    if(nb_arg<argc-1){
      strcat(cmdLine," ");
    }
  }
  //printf("\nCommand line: '%s'\n",cmdLine);

  // copy argv values to the respective global arguments
  for(nb_arg=1; nb_arg<argc; nb_arg++){

    if(strcmp(argv[nb_arg],"-h")==0){
      printf("\n%s",usage);
      return(24);
    }

    if(strcmp(argv[nb_arg],"-p1")==0){
      nrg_file1=argv[nb_arg+1];
    }

    if(strcmp(argv[nb_arg],"-p2")==0){
      nrg_file2=argv[nb_arg+1];
    }

    if(strcmp(argv[nb_arg],"-pp")==0){
      pairwiseF=argv[nb_arg+1];
    }

    if(strcmp(argv[nb_arg],"-w")==0){
      wrfn=1;
    }

    if(strcmp(argv[nb_arg],"-q")==0){
      char* pch = strtok (argv[nb_arg+1],",");
      while (pch != NULL){
        list1.push_back(atoi(pch));
        pch = strtok (NULL, ",");
      }

      pch = strtok (argv[nb_arg+2],",");
      while (pch != NULL){
        list2.push_back(atoi(pch));
        pch = strtok (NULL, ",");
      }

      if(list1.size() != list2.size()){
        cout<<"Lists of atoms to superimpose must have same nb of atoms."<<endl;
        return(24);
      }
    }

    if(strcmp(argv[nb_arg],"-d")==0){
      sscanf(argv[nb_arg+1],"%f",&dDist);
    }

    if(strcmp(argv[nb_arg],"-dca")==0){
      sscanf(argv[nb_arg+1],"%f",&ca_dDist);
    }

    if(strcmp(argv[nb_arg],"-dps")==0){
      sscanf(argv[nb_arg+1],"%f",&ps_dDist);
    }

    if(strcmp(argv[nb_arg],"-dcn")==0){
      sscanf(argv[nb_arg+1],"%f",&neibr_dDist);
    }

    if(strcmp(argv[nb_arg],"-s")==0){
      if(strcmp(argv[nb_arg+1],"1")==0){
        bkAll=1;
      }else if(strcmp(argv[nb_arg+1],"0")==0){
        bkAll=0;
      }
    }

    if(strcmp(argv[nb_arg],"-n")==0){
      commonInt=atoi(argv[nb_arg+1]);
    }

    if(strcmp(argv[nb_arg],"-e")==0){
      emptOut=1;
    }

    if(strcmp(argv[nb_arg],"-a")==0){
      maxCliques=atoi(argv[nb_arg+1]);
    }

    if(strcmp(argv[nb_arg],"-j")==0){
      jttt=atoi(argv[nb_arg+1]);
    }

    if(strcmp(argv[nb_arg],"-wc")==0){
      wc=1;
    }

    if(strcmp(argv[nb_arg],"-m")==0){
      maxNodes=atoi(argv[nb_arg+1]);
    }

    if(strcmp(argv[nb_arg],"-c")==0){
      char* pch = strtok (argv[nb_arg+1],",");
      while (pch != NULL){
        steps.push_back(atoi(pch));
        if(atoi(pch)>=0 && cg_start<0){ cg_start=atoi(pch); }
        pch = strtok (NULL, ",");
      }
    }

    if(strcmp(argv[nb_arg],"-o")==0){
      sprintf(outbase,"%s",argv[nb_arg+1]);
    }

    if(strcmp(argv[nb_arg],"-l1")==0){
      rnc1=argv[nb_arg+1];
    }

    if(strcmp(argv[nb_arg],"-l2")==0){
      rnc2=argv[nb_arg+1];
    }

    if(strcmp(argv[nb_arg],"-pr")==0){
      printDetails=1;
    }
  }

  return(0);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/

//This Function opens the two mif file and gets the two prefix of the proteins
//and also looks if they have the same number of probes
int get_info(string str1, string str2){
  int i;
  const char* prefix1; //Mif file 1 prefix
  const char* prefix2; //Mif file 2 prefix
  const char* prefix3; //Pairwise prefix
  const char* p1_pre; //Protein file 1 prefix
  const char* p2_pre; //Protein file 2 prefix

  //Store Mif 1 filename without the .isomif in prefix1
  string empt ("");
  string pre1;
  size_t found1;
  found1=str1.find(".mif");
  if(found1!=string::npos){
    for(i=found1-1; i>=0; i--){
      if(str1.at(i)=='/'){
        pre1=str1.substr(i+1,found1-i-1);
        break;
      }
    }
  }
  if(pre1.compare(empt)==0){
    pre1=str1.substr(0,found1-i-1);
  }
  prefix1 = pre1.c_str();
  tag1=pre1;

  //Store Mif 2 filename without the .isomif in prefix1
  string pre2;
  size_t found2;
  found2=str2.find(".mif");
  if(found2!=string::npos){
    for(i=found2-1; i>=0; i--){
      if(str2.at(i)=='/'){
       pre2=str2.substr(i+1,found2-i-1);
       break;
     }
   }
  }
  if(pre2.compare(empt)==0){
    pre2=str2.substr(0,found2-i-1);
  }
  prefix2 = pre2.c_str();
  tag2=pre2;

  //Create output file
  if(!strcmp(outbase,"")){ //If no outbase name is defined, create one
    // sprintf(out_file,"./%s_match_%s.pdb",prefix1,prefix2);
    sprintf(out_file,"./%s_match_%s",prefix1,prefix2);
  }else{
    // sprintf(out_file,"%s%s_match_%s.pdb",outbase,prefix1,prefix2);
    sprintf(out_file,"%s%s_match_%s",outbase,prefix1,prefix2);
  }


  return(0);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void rem_spaces(char *ch){
    char ch2[strlen(ch)+1];
    int i,j=0;
    for(i=0;i<=strlen(ch);i++){
      if(ch[i]!=' '){
        ch2[j] = ch[i];
        j++;
      }
    }
    strcpy(ch,ch2);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/

long long ConnID(int id1, int id2){

  long long id=0;
  long long x = (long long)id1;
  long long y = (long long)id2;
  long long nn = (long long)numNodes;

  if(id1 < id2){
    id = (( (x*nn) + (y + 1) ) - ( ((x + 1) * x) / 2 ))-1;
  }else{
    id = (( (y*nn) + (x + 1) ) - ( ((y + 1) * y) / 2 ))-1;
  }
  return(id);
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void setJTT(int threshold){
  int i,j;
  int jtt_ranks[20][20] =
    {
      {0,8,10,7,14,10,6,4,15,12,9,11,13,15,5,2,1,17,16,3},
      {6,0,10,14,9,2,11,3,4,12,8,1,13,16,7,5,7,11,15,12},
      {7,11,0,2,15,9,8,6,5,10,13,3,15,16,14,1,4,17,12,13},
      {4,10,2,0,13,8,1,3,6,11,11,9,12,13,10,5,7,14,9,8},
      {4,3,10,13,0,14,15,4,9,11,8,14,12,6,12,1,7,9,2,5},
      {7,3,8,10,17,0,1,11,4,15,6,2,13,17,5,8,9,16,14,12},
      {5,10,7,1,15,2,0,4,12,12,11,3,13,14,11,8,9,15,14,6},
      {2,4,6,5,11,12,3,0,15,16,14,9,16,17,10,1,8,13,16,7},
      {10,3,4,6,13,1,12,11,0,14,7,9,15,12,5,8,9,16,2,15},
      {6,10,8,13,16,15,12,15,15,0,2,9,4,5,14,7,3,17,11,1},
      {8,9,15,16,15,7,14,15,11,1,0,12,4,3,5,6,10,13,12,2},
      {6,1,3,9,15,2,3,7,11,11,10,0,8,14,11,5,4,15,13,12},
      {5,7,11,13,15,9,12,12,14,2,1,6,0,10,13,8,4,16,15,3},
      {7,13,12,14,6,14,13,12,9,4,1,14,8,0,9,3,11,10,2,5},
      {2,6,13,12,16,5,11,8,7,14,3,10,15,13,0,1,4,17,16,9},
      {1,7,3,10,8,13,16,4,17,14,6,11,18,9,5,0,2,19,15,12},
      {1,8,4,11,13,10,11,9,12,3,9,6,7,15,6,2,0,16,14,5},
      {9,1,12,13,5,8,10,3,13,11,2,10,10,6,12,4,9,0,5,7},
      {12,11,5,6,4,11,15,13,2,8,7,14,15,1,14,3,9,11,0,10},
      {2,12,14,10,11,13,7,6,17,1,3,13,4,9,11,8,5,16,15,0}
    };

  for(i=0;i<20;i++){
    for(j=0;j<20;j++){
      jtt[i][j] =0;
      if(jtt_ranks[i][j] < threshold) jtt[i][j] = 1;
    }
  }
  return;
}
/***********************************************************************/
/*        1         2         3         4         5         6         7*/
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int nam2n(string rnam){
  if(rnam.compare("ALA") == 0) return 0;
  if(rnam.compare("ARG") == 0) return 1;
  if(rnam.compare("ASN") == 0) return 2;
  if(rnam.compare("ASP") == 0) return 3;
  if(rnam.compare("CYS") == 0) return 4;
  if(rnam.compare("GLN") == 0) return 5;
  if(rnam.compare("GLU") == 0) return 6;
  if(rnam.compare("GLY") == 0) return 7;
  if(rnam.compare("HIS") == 0) return 8;
  if(rnam.compare("ILE") == 0) return 9;
  if(rnam.compare("LEU") == 0) return 10;
  if(rnam.compare("LYS") == 0) return 11;
  if(rnam.compare("MET") == 0) return 12;
  if(rnam.compare("PHE") == 0) return 13;
  if(rnam.compare("PRO") == 0) return 14;
  if(rnam.compare("SER") == 0) return 15;
  if(rnam.compare("THR") == 0) return 16;
  if(rnam.compare("TRP") == 0) return 17;
  if(rnam.compare("TYR") == 0) return 18;
  if(rnam.compare("VAL") == 0) return 19;

  return -1;

}

// ##### FAKE GRAPH ####### //
    // numNodes=23;
    // vector<node> graph2;
    // for(int a=0; a<numNodes; a++){
    //   node newNode;
    //   newNode.id = a;
    //   newNode.uid = a;
    //   graph2.push_back(newNode);
    // }
    // int nvsize = (((numNodes*numNodes)-numNodes)/2)+numNodes;
    // bool* conn2=new bool[nvsize];
    // for(i=0; i<graph2.size(); ++i){
    //   conn2[ConnID(i,i)]=1;
    //   for(j=0; j<i; ++j){
    //     conn2[ConnID(i,j)]=0;
    //     // cout<<i<<" "<<j<<" "<< ConnID(i,j)<<endl;
    //   }
    // }
    // conn2[ConnID(1,7)]=1;
    // conn2[ConnID(1,3)]=1;
    // conn2[ConnID(1,6)]=1;
    // conn2[ConnID(1,5)]=1;
    // conn2[ConnID(3,5)]=1;
    // conn2[ConnID(3,6)]=1;
    // conn2[ConnID(6,5)]=1;
    // conn2[ConnID(5,7)]=1;
    // conn2[ConnID(0,2)]=1;
    // conn2[ConnID(0,4)]=1;
    // conn2[ConnID(2,4)]=1;
    // conn2[ConnID(0,3)]=1;
    // conn2[ConnID(0,6)]=1;
    // conn2[ConnID(0,7)]=1;    
    // conn2[ConnID(6,11)]=1;
    // conn2[ConnID(6,10)]=1;
    // conn2[ConnID(6,8)]=1;
    // conn2[ConnID(10,8)]=1; 
    // conn2[ConnID(11,8)]=1;
    // conn2[ConnID(11,10)]=1;  
    // conn2[ConnID(11,9)]=1;
    // conn2[ConnID(9,10)]=1;  
    // conn2[ConnID(9,8)]=1;  
    // conn2[ConnID(9,6)]=1;  
    // conn2[ConnID(12,8)]=1;  
    // conn2[ConnID(12,6)]=1;
    // conn2[ConnID(12,11)]=1;  
    // conn2[ConnID(12,9)]=1;  
    // conn2[ConnID(12,10)]=1;  
    // conn2[ConnID(7,13)]=1;
    // conn2[ConnID(7,14)]=1;  
    // conn2[ConnID(7,15)]=1;  
    // conn2[ConnID(7,16)]=1;    
    // conn2[ConnID(13,14)]=1;  
    // conn2[ConnID(13,15)]=1;  
    // conn2[ConnID(13,16)]=1; 
    // conn2[ConnID(16,14)]=1;  
    // conn2[ConnID(16,15)]=1;
    // conn2[ConnID(14,15)]=1; 
    // conn2[ConnID(17,7)]=1;  
    // conn2[ConnID(17,13)]=1; 
    // conn2[ConnID(17,16)]=1;  
    // conn2[ConnID(17,15)]=1;
    // conn2[ConnID(17,14)]=1; 
    // conn2[ConnID(7,18)]=1;
    // conn2[ConnID(7,19)]=1;
    // conn2[ConnID(20,8)]=1;
    // conn2[ConnID(20,9)]=1;
    // conn2[ConnID(20,6)]=1;
    // conn2[ConnID(20,10)]=1; 

    // cout<<"Entering Bron Kerbosch"<<endl;
    // bk(cg,graph2,conn2);