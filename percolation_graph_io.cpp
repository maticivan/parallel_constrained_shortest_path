

//*******************************************************************************************************
//******************                                                                                *****
//******************  A parallel algorithm for constrained shortest path problem                    *****
//******************                                                                                *****
//******************  Created by Ivan Matic                                                         *****
//******************  http://www.imomath.com/maticivan                                              *****
//******************                                                                                *****
//******************  License: GNU General Public License, version 3 (GPL-3.0)                      *****
//******************                                                                                *****
//******************  The Software is provided "as is," with all faults, defects and errors,        *****
//******************  and without warranty of any kind.                                             *****
//******************  Licensor does not warrant that the Software will be free of bugs, errors,     *****
//******************  viruses or other defects, and Licensor shall have no liability of any         *****
//******************  kind for the use of or inability to use the software, the software content    *****
//******************  or any associated service.                                                    *****
//*******************************************************************************************************

#ifndef PERCOLATION_IO
#define PERCOLATION_IO

#define myint int
#define mysint int



#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <vector>
#include <limits>
#include <ctime>
#include <iomanip>
#include <math.h>
#include <functional>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

using namespace std;
using namespace boost;



struct EdgeProperties {
    myint eTime, eCost;
    mysint eState;// 0-passive; 1-used; 2-just used; 3-active; 4-active and activated in this cycle
    myint eSource; // vName of the source
    myint edgeId1;
    myint edgeId2;
};

struct VertexProperties {
    myint vName;
    myint vCost;
    mysint belongsToD;
    mysint belongsToC;
    myint originalName;
    mysint vState;// 0- inactive; 1-active
};


typedef adjacency_list < vecS, vecS, undirectedS, VertexProperties, EdgeProperties> Graph;
typedef graph_traits<Graph>::vertex_descriptor VertexDesc;
typedef graph_traits<Graph>::edge_descriptor EdgeDesc;
typedef graph_traits<Graph>::edge_iterator EdgeIt;
typedef graph_traits<Graph>::adjacency_iterator AdjIt;

typedef graph_traits<Graph>::vertex_iterator VtxIt;



typedef std::numeric_limits< double > dbl;




int getNextInteger(char * memblock, streampos size, streampos *pos){
    
    int current=0;
    int mult=1;
    
    while ((*pos<size) &&(memblock[*pos]!='-')&&(( memblock[*pos]<'0')||(memblock[*pos]>'9'))){
        *pos= (*pos)+((streampos)1);
    }
    if(memblock[*pos]=='-'){
        mult=-1;
        *pos= (*pos)+((streampos)1);
    }
    while ((*pos<size) &&(( memblock[*pos]>='0')&&(memblock[*pos]<='9'))){
        current= 10 * current;
        current+= (int)(memblock[*pos]-'0');
        *pos=*pos+((streampos)1);
    }
    current=current*mult;
    return current;
}
myint randomNumber(myint randNk){
    myint random_integer=rand()%randNk;
    return random_integer;
}


myint getSetupFromFile(string filename, myint *numThreads, myint *dimensionOld, myint *dimensionNew, myint *maxTime, myint * maxCost, myint *n , myint *boundaryVertCost, myint *boundaryChoice, myint *debuggingMode, myint *cpuForce){
    
    streampos size; streampos *position;
    position=new streampos;
    *position=0;
    char * memblock;
    ifstream ifile(filename,ios::in|ios::binary|ios::ate);
    myint modulo;
    if (ifile.is_open())
    {
        size = ifile.tellg();
        memblock = new char [size];
        ifile.seekg (0, ios::beg);
        ifile.read (memblock, size);
        ifile.close();
        *numThreads=getNextInteger(memblock,size,position);
        *dimensionOld=getNextInteger(memblock,size,position);
        *dimensionNew=getNextInteger(memblock,size,position);
        *maxTime=getNextInteger(memblock,size,position);
        *maxCost=getNextInteger(memblock,size,position);
        *n=getNextInteger(memblock,size,position);
        *boundaryVertCost=getNextInteger(memblock,size,position);
        *boundaryChoice=getNextInteger(memblock,size,position);
        
        *debuggingMode=getNextInteger(memblock,size,position);
        *cpuForce=getNextInteger(memblock,size,position);
        delete[] memblock;
        
    }
    return 1;
}




int readFromFile(string filename, Graph **pMyGraph){
    
    streampos size; streampos *position;
    position=new streampos;
    *position=0;
    Graph *myGraph;
    
    char * memblock;
    ifstream ifile(filename,ios::in|ios::binary|ios::ate);
    EdgeDesc edge1;
    int counter, int1,int2,int3,int4,tempSource,tempState, tempBelongsToC;
    int tempBelongsToD, tempOriginalName, minus1;
    minus1=-1;
    if (ifile.is_open())
    {
        size = ifile.tellg();
        memblock = new char [size];
        ifile.seekg (0, ios::beg);
        ifile.read (memblock, size);
        ifile.close();
        
        int fileType;
        fileType=getNextInteger(memblock,size,position);
        
        int1=getNextInteger(memblock, size, position);
        
        myGraph=new Graph(int1);
        *pMyGraph=myGraph;
        VtxIt begV,endV;
        AdjIt begE,endE;
        int1=getNextInteger(memblock,size,position);
        if(int1==-7){
            int1 = getNextInteger(memblock,size,position);// number of edges
            int1 = getNextInteger(memblock,size,position);
        }
        tie(begV,endV)=vertices(*myGraph);
        counter=0;
        while((int1!=-1)&&(begV!=endV)){
            put(&VertexProperties::vName,*myGraph,*begV,counter);

            put(&VertexProperties::vCost,*myGraph,*begV,int1);
            put(&VertexProperties::originalName,*myGraph,*begV,minus1);
            tempState=0;
            tempBelongsToD=0;tempBelongsToC=0;tempOriginalName=-1;
            if(int1>0){
                tempState=1;
                tempBelongsToD=1;
            }

            if(fileType==2){
                int1=getNextInteger(memblock,size,position);//Name - already know it
                tempState=getNextInteger(memblock,size,position);//Cost - already know it
                tempState=getNextInteger(memblock,size,position);
                tempBelongsToD= getNextInteger(memblock,size,position);
                tempBelongsToC=getNextInteger(memblock,size,position);
                tempOriginalName= getNextInteger(memblock,size,position);
            }
            put(&VertexProperties::vState,*myGraph,*begV,tempState);
            
            put(&VertexProperties::belongsToD,*myGraph,*begV,tempBelongsToD);
            put(&VertexProperties::belongsToC,*myGraph,*begV,tempBelongsToC);
            

            
            int1=getNextInteger(memblock,size,position);
            
            counter++;begV++;
            
        }
        while(begV!=endV){
            int1=0;
            put(&VertexProperties::vName,*myGraph,*begV,counter);
            put(&VertexProperties::vCost,*myGraph,*begV,int1);
            put(&VertexProperties::originalName,*myGraph,*begV,minus1);
            put(&VertexProperties::belongsToC,*myGraph,*begV,0);
            put(&VertexProperties::belongsToD,*myGraph,*begV,0);
            put(&VertexProperties::vState,*myGraph,*begV,0);
            counter++;begV++;
        }

        
        int1=getNextInteger(memblock,size,position);
        if(int1!=-1){
            int2=getNextInteger(memblock,size,position);
            int3=getNextInteger(memblock,size,position);
            int4=getNextInteger(memblock,size,position);
        }
        while(int1!=-1){
            if(!(edge(int1,int2,*myGraph).second)){
                
                edge1=add_edge(int1,int2,*myGraph).first;
                put(&EdgeProperties::eTime,*myGraph,edge1,int3);
                put(&EdgeProperties::eCost,*myGraph,edge1,int4);
                tempState=0;tempSource=-1;
                if(get(&VertexProperties::vCost,*myGraph,int1)>0){
                    tempState=3;tempSource=int1;
                }
                if(get(&VertexProperties::vCost,*myGraph,int2)>0){
                    tempState=3;tempSource=int2;
                }
                if (fileType==2){
                    tempState=getNextInteger(memblock,size,position);
                    tempSource=getNextInteger(memblock,size,position);
                }
                put(&EdgeProperties::eState,*myGraph,edge1,tempState);
                put(&EdgeProperties::eSource,*myGraph,edge1,tempSource);
            }
            else{
                if(fileType==2){
                    tempState=getNextInteger(memblock,size,position);
                    tempSource=getNextInteger(memblock,size,position);
                }
            }
            int1=getNextInteger(memblock,size,position);
            if(int1!=-1){
                int2=getNextInteger(memblock,size,position);
                int3=getNextInteger(memblock,size,position);
                int4=getNextInteger(memblock,size,position);

            }
        }
        
        delete[] memblock;
    }
    
    
    
    return 1;
}

//the last parameter Names tells the kind of output we want.
//If names==0, then the output will be of the form that can used as an input of a graph in some future
//execution of the program. Such an output will be minimalistic and will not have good display quality.
//If names==1, then the output will be of display quality for debugging and checking the program.

int printToFile(string filename, Graph *myGraph, int names){
    ofstream mfile;
    mfile.open(filename);
    
    VtxIt begV,endV;
    AdjIt begE,endE;
    EdgeDesc edge1;
    int number2=2;
    mfile<<number2<<endl;
    mfile<<num_vertices(*myGraph)<<endl;
    mfile<< "-7 "<<num_edges(*myGraph)<<endl;
    tie(begV,endV)=vertices(*myGraph);
    for(tie(begV,endV)=vertices(*myGraph);begV!=endV;++begV){
        mfile<<get(&VertexProperties::vCost,*myGraph,*begV);
        if(names==1){
            mfile<<" Name: "<<get(&VertexProperties::vName,*myGraph,*begV)<<" ";
            mfile<<" Cost: "<<get(&VertexProperties::vCost,*myGraph,*begV)<<" ";
            mfile<<" State: "<<get(&VertexProperties::vState,*myGraph,*begV)<<" ";
            mfile<<" In D: "<<get(&VertexProperties::belongsToD,*myGraph,*begV)<<" ";
            mfile<<" In C: "<<get(&VertexProperties::belongsToC,*myGraph,*begV)<<" ";
            mfile<<" originalName: "<<get(&VertexProperties::originalName,*myGraph,*begV)<<" ";

        }
        mfile<<endl;
    }
    int minus1=-1;
    mfile<<minus1<<endl;
    for(tie(begV,endV)=vertices(*myGraph);begV!=endV;++begV){
        tie(begE,endE)=adjacent_vertices(*begV,*myGraph);
        if((begE==endE)&&(names==1)){
            mfile<<*begV<<" ";
            if(names==1){
                mfile<<get(&VertexProperties::vName,*myGraph,*begV)<<" ";
            }
            mfile<<"no connections"<<endl;
        }
        for(tie(begE,endE)=adjacent_vertices(*begV,*myGraph);begE!=endE;++begE){
            mfile<<"("<<*begV<<","<<*begE<<") Time: ";
            edge1=edge(*begV,*begE,*myGraph).first;
            mfile<<get(&EdgeProperties::eTime,*myGraph,edge1)<<" Cost: ";
            mfile<<get(&EdgeProperties::eCost,*myGraph,edge1);

            mfile<<" State: "<<get(&EdgeProperties::eState,*myGraph,edge1);
            mfile<<" Source: "<<get(&EdgeProperties::eSource,*myGraph,edge1);
            mfile<<endl;
        }
    }
    
    mfile<<minus1<<" "<<minus1<<" "<<minus1<<" "<<minus1<<" "<<minus1<<" "<<minus1<<endl;
    
    mfile.close();
    return 1;
}

myint vertexName(myint *coordinates, myint dimension, myint n){
    myint tempExp=0;
    myint tempResult=0;
    while(tempExp<dimension){
        tempExp++;
        tempResult=tempResult * n+ coordinates[dimension-tempExp];
    }
    return tempResult;
}
myint incrementVertex(myint *coordinates, myint dimension, myint n){
    myint i=0, stillGoing=1;
    while(stillGoing==1){
        if(coordinates[i]<n-1){
            coordinates[i]=coordinates[i]+1;
            stillGoing=0;
        }
        else{
            if(i<dimension){
                coordinates[i]=0;
                i++;
            }
            else{
                stillGoing=2;
            }
        }
    }
    return stillGoing;
}

myint createBigGraph(Graph **pMyGraph, myint d, myint n, myint maxTime, myint maxCost, myint boundaryVertCost, myint boundaryChoice){
    
    //myint boundaryVertCost=19;
    myint boundaryVertCostTemp=0;
    Graph *myGraph;
    

    EdgeDesc edge1;
    myint tempSource,tempState,minus1;
    minus1=-1;
    myint nToTheD=1;
    for(myint i=0;i<d;i++){
        nToTheD*=n;
    }
    myGraph=new Graph(nToTheD);
    *pMyGraph=myGraph;
    VtxIt begV,endV;
    AdjIt begE,endE;
    tie(begV,endV)=vertices(*myGraph);
    myint *x; myint *y;
    x=new myint[d];
    y=new myint[d];
    for (myint i=0;i<d;i++){
        x[i]=0;
    }
    myint boundaryIndicator, centerIndicator;
    myint nameVariable, nameVariableY, int3, int4;
 
    while(begV!=endV){
        nameVariable=vertexName(x,d,n);
        put(&VertexProperties::vName,*myGraph,*begV,nameVariable);
        boundaryIndicator=0;boundaryVertCostTemp=0; centerIndicator=1;
        
        for(myint i=0;i<d;i++){
            if((x[i]==0)||(x[i]==n-1)){
                boundaryIndicator=1;
                boundaryVertCostTemp=boundaryVertCost;
            }
            if( (2*x[i]>n+1)||(2*x[i]<n-1)){
                centerIndicator=0;
            }
        }
        if(boundaryChoice==1){
            for(myint i=0;i<d;i++){
                if(x[i]!=0){
                    boundaryIndicator=0;
                    boundaryVertCostTemp=0;
                }
            }
        }
        put(&VertexProperties::vState,*myGraph,*begV,boundaryIndicator);
        put(&VertexProperties::vCost,*myGraph,*begV,boundaryVertCostTemp);
        put(&VertexProperties::originalName,*myGraph,*begV,minus1);
        if(centerIndicator==1){
            boundaryIndicator=2;
        }
        put(&VertexProperties::belongsToD,*myGraph,*begV,boundaryIndicator);
        put(&VertexProperties::belongsToC,*myGraph,*begV,0);
        
        begV++; incrementVertex(x,d,n);
    }
    
    tie(begV,endV)=vertices(*myGraph);
    while(begV!=endV){
        nameVariable=vertexName(x,d,n);
        
        for(myint i=0;i<d;i++){
            if(x[i]<n-1){
                for(myint j=0;j<d;j++){
                    y[j]=x[j];
                }
                y[i]=y[i]+1;
                nameVariableY=vertexName(y,d,n);
                if(!(edge(nameVariable ,nameVariableY,*myGraph).second)){
                    
                    edge1=add_edge(nameVariable,nameVariableY,*myGraph).first;
                    int3=randomNumber(maxTime-1)+1;
                    int4=randomNumber(maxCost-1)+1;
                    put(&EdgeProperties::eTime,*myGraph,edge1,int3);
                    put(&EdgeProperties::eCost,*myGraph,edge1,int4);
                    
                    
                    
                    tempState=0;tempSource=-1;
                    if(get(&VertexProperties::vCost,*myGraph,nameVariable)>0){
                        tempState=3;tempSource=nameVariable;
                    }
                    if(get(&VertexProperties::vCost,*myGraph,nameVariableY)>0){
                        tempState=3;tempSource=nameVariableY;
                    }
                    put(&EdgeProperties::eState,*myGraph,edge1,tempState);
                    put(&EdgeProperties::eSource,*myGraph,edge1,tempSource);
                }
            }
        }
        
        begV++; incrementVertex(x,d,n);
    }
    delete[] x; delete[] y;
 
    
    return 1;
}
#endif
