
//*******************************************************************************************************
//******************                                                                                *****
//******************  A parallel algorithm for constrained shortest path problem                    *****
//******************                                                                                *****
//******************  Created by Ivan Matic                                                         *****
//******************  License: GNU General Public License, version 3 (GPL-3.0)                      *****
//******************                                                                                *****
//******************  The Software is provided "as is," with all faults, defects and errors,        *****
//******************  and without warranty of any kind.                                             *****
//******************  Licensor does not warrant that the Software will be free of bugs, errors,     *****
//******************  viruses or other defects, and Licensor shall have no liability of any         *****
//******************  kind for the use of or inability to use the software, the software content    *****
//******************  or any associated service.                                                    *****
//*******************************************************************************************************


// compile with C++ -I /usr/local/boost_1_58_0 -framework opencl percolation.cpp
// compile with C++ -I /usr/local/boost_1_59_0 -framework opencl percolation.cpp


#include "percolation_graph_io.cpp"

#include "percolation_parallel.cpp"



int main(){
    
    srand((unsigned)time(0));
    
    
    myint oldDim, dim,n, maxTime, maxCost, numTh, debuggingMode, cpuForce, threadForce, avoidBoost, bVC, bChoice, tFP1, tFP2;
    getSetupFromFile("setup.txt", &numTh, &oldDim, &dim, &maxTime, &maxCost, &n , &bVC, &bChoice, &debuggingMode, &cpuForce);
    
    Graph **pgr,**pgr2;
    pgr=new Graph*;
    pgr2=new Graph*;
    int totalTime=0;
    int glj;
 
        readFromFile("example.txt",pgr);
        copyGraph(*pgr,pgr2);
    
    shortestPathParallel(*pgr2);
    
 
    delete *pgr;
    delete pgr;
    delete *pgr2;
    delete pgr2;
    
 
    pgr=new Graph*;
    
    
    
    
    createBigGraph(pgr, dim,n,maxTime, maxCost, bVC, bChoice);
    printToFile("newRandomGraph.txt",*pgr,1);
    delete *pgr;
    delete pgr;
    
    
    
    
    
    
    
    
    return 1;
    
}
