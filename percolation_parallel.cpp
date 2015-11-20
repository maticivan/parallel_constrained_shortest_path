

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


mysint rnk(mysint randNk){
    mysint random_integer=rand()%randNk;
    return random_integer;
}
int copyGraph (Graph *sourceG, Graph **destination){
    
    Graph *myGraph;
    
    
    myGraph =new Graph(*sourceG);
    *destination=myGraph;
    
    return 1;
}





cl_context CreateContext(myint forceCPU)
{
    cl_int errNum;
    cl_uint numPlatforms;
    cl_platform_id firstPlatformId;
    cl_context context = NULL;
    
    
    errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
    if (errNum != CL_SUCCESS || numPlatforms <= 0)
    {
        std::cerr << "No OpenCL platforms." << std::endl;
        return NULL;
    }

    cl_context_properties contextProperties[] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)firstPlatformId,
        0
    };
    
    //Select GPU:
    
 
    if(forceCPU==0){
        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU,NULL, NULL, &errNum);
        if (errNum != CL_SUCCESS)
        {
            std::cout << "No GPU context. Trying CPU." << std::endl;
            context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
                                              NULL, NULL, &errNum);
            if (errNum != CL_SUCCESS)
            {
                std::cerr << "Failed to create an OpenCL GPU or CPU context." << std::endl;
                return NULL;
            }
        }
    }
    //Comparison: Select CPU:

    if(forceCPU==1){
        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
                                          NULL, NULL, &errNum);
        if (errNum != CL_SUCCESS)
        {
            std::cerr << "Failed to create an OpenCL GPU or CPU context." << std::endl;
            return NULL;
        }
    }
    
    

    return context;
}


cl_command_queue CreateCommandQueue(cl_context context, cl_device_id *device)
{
    cl_int errNum;
    cl_device_id *devices;
    cl_command_queue commandQueue = NULL;
    size_t deviceBufferSize = -1;
    
 
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceBufferSize);

    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Failed call to clGetContextInfo(...,GL_CONTEXT_DEVICES,...)";
        return NULL;
    }
    
    if (deviceBufferSize <= 0)
    {
        std::cerr << "No devices available.";
        return NULL;
    }
 
    myint numDevices=(myint)(deviceBufferSize / sizeof(cl_device_id));
    devices = new cl_device_id[deviceBufferSize / sizeof(cl_device_id)];
    
 
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceBufferSize, devices, NULL);
    if (errNum != CL_SUCCESS)
    {
        delete [] devices;
        std::cerr << "Failed to get device IDs";
        return NULL;
    }
    
    // If I am at computer at work then deviceToChoose should be 0; if I am at home it should be 1;
    
    myint deviceToChoose=numDevices-1;
   // deviceToChoose=0;

    
    commandQueue = clCreateCommandQueue(context, devices[deviceToChoose], 0, NULL);

    if (commandQueue == NULL)
    {
        delete [] devices;
        std::cerr << "Failed to create commandQueue for device 0";
        return NULL;
    }
    
    *device = devices[deviceToChoose];
    delete [] devices;
    return commandQueue;
}


cl_program CreateProgram(cl_context context, cl_device_id device, const char* fileName)
{
    cl_int errNum;
    cl_program program;
    
    std::ifstream kernelFile(fileName, std::ios::in);
    if (!kernelFile.is_open())
    {
        std::cerr << "Failed to open file for reading: " << fileName << std::endl;
        return NULL;
    }
    
    std::ostringstream oss;
    oss << kernelFile.rdbuf();
    
    std::string srcStdStr = oss.str();
    const char *srcStr = srcStdStr.c_str();
    program = clCreateProgramWithSource(context, 1,
                                        (const char**)&srcStr,
                                        NULL, NULL);
    if (program == NULL)
    {
        std::cerr << "Failed to create CL program from source." << std::endl;
        return NULL;
    }
    
    errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (errNum != CL_SUCCESS)
    {
 
        char buildLog[16384];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                              sizeof(buildLog), buildLog, NULL);
        
        std::cerr << "Error in kernel: " << std::endl;
        std::cerr << buildLog;
        clReleaseProgram(program);
        return NULL;
    }
    
    return program;
}

//*************************************************************************************
//***************************      mergeSortDec
//*************************************************************************************
//
//
// The function performs these two steps:
// 1) Sorts the sequence sequenceToBeSorted of length lengthOfTheSequenceToBeSorted in decreasing order
// 2) The obtained sequence is then merged with alreadySortedSequence of length lengthOfSortedSequence
//    that is assumed to be already sorted in decreasing order.
// The merged sequence will be stored in sequenceToBeSorted and its final length is lengthOfTheSequenceToBeSorted+lengthOfSortedSequence.


myint mergeSortDec(myint *sequenceToBeSorted, myint lengthOfTheSequenceToBeSorted,  myint *alreadySortedSequence=NULL, myint lengthOfSortedSequence=0){
    
    if(lengthOfTheSequenceToBeSorted>1){
        myint middle=(lengthOfTheSequenceToBeSorted+1)/2;
        mergeSortDec(sequenceToBeSorted+middle,lengthOfTheSequenceToBeSorted-middle);
        mergeSortDec(sequenceToBeSorted,middle,sequenceToBeSorted+middle,lengthOfTheSequenceToBeSorted-middle);
    }
    if(lengthOfSortedSequence>0){
        myint length=lengthOfTheSequenceToBeSorted+lengthOfSortedSequence;
        myint* help;
        help=new myint[length];
        myint r=0, rL=0,rR=0;
        while(r<length){
            if(rL>=lengthOfTheSequenceToBeSorted){
                help[r]=alreadySortedSequence[rR];
                rR++;
            }
            else{
                if(rR>=lengthOfSortedSequence){
                    help[r]=sequenceToBeSorted[rL];
                    rL++;
                }
                else{
                    if(sequenceToBeSorted[rL]>alreadySortedSequence[rR]){
                        help[r]=sequenceToBeSorted[rL];
                        rL++;
                    }
                    else{
                        help[r]=alreadySortedSequence[rR];
                        rR++;
                    }
                }
            }
            r++;
        }
        for(r=0;r<length;r++){
            sequenceToBeSorted[r]=help[r];
        }
        delete[] help;
    }
    return 0;
}

//*************************************************************************************
//***************************      mergeSortDecBlock
//*************************************************************************************
//
//
// mergeSortDecBlock is used for sorting the sequence of phantoms.
// It sorts the sequence unsortedS in decreasing order and merges it with the sequence alreadySorted
// It differs from mergeSortDecBlock in that the sequences are blocks of size blockSize and the relation ">" considers the first 5 terms and
// is defined in the following way:
// (a_0,a_1,a_2,a_3,a_4)>(b_0,b_1,b_2,b_3,b_4) if one of the following 3 cases occurs:
// Case 1: a_1>b_1
// Case 2: a_1=b_1 and a_3>b_3
// Case 3: a_1=b_1 and a_3=b_3 and (a_2-a_4)>(b_2-b_4)
// If neither a>b nor b>a then a and b are called equivalent and eventually one of them will be erased from the sequence.

myint mergeSortDecBlock(myint *unsortedS, myint lengthOfUnsorted, myint blockSize,  myint *alreadySorted=NULL, myint lengthOfSorted=0){
    myint counter,comparisonInd;
    if(lengthOfUnsorted>1){
        myint middle=(lengthOfUnsorted+1)/2;
        mergeSortDecBlock(unsortedS+(middle*blockSize),lengthOfUnsorted-middle,blockSize);
        mergeSortDecBlock(unsortedS,middle,blockSize,unsortedS+(middle*blockSize),lengthOfUnsorted-middle);
    }
    if(lengthOfSorted>0){
        myint length=lengthOfUnsorted+lengthOfSorted;
        myint* help;
        help=new myint[length * blockSize];
        myint r=0, rL=0,rR=0;
        while(r<length){
            if(rL>=lengthOfUnsorted){
                for(counter=0;counter<blockSize;counter++){
                    help[r*blockSize+counter]=alreadySorted[rR*blockSize+counter];
                }
                rR++;
            }
            else{
                if(rR>=lengthOfSorted){
                    for(counter=0;counter<blockSize;counter++){
                        help[r*blockSize+counter]=unsortedS[rL*blockSize+counter];
                    }
                    rL++;
                }
                else{
                    comparisonInd=0;
                    //comparisonInd should be =1 if unsortedS[rL]>alreadySorted[rR]
                    //                        0 otherwise
                    
                    if(unsortedS[rL*blockSize+1]>alreadySorted[rR*blockSize+1]){
                        comparisonInd=1;
                    }
                    else{
                        if(unsortedS[rL * blockSize+1]==alreadySorted[rR*blockSize+1]){
                            if(unsortedS[rL * blockSize+3]>alreadySorted[rR*blockSize+3]){
                                comparisonInd=1;
                            }
                            else{
                                if(unsortedS[rL * blockSize+3]==alreadySorted[rR*blockSize+3]){
                                    if(unsortedS[rL * blockSize+2]-unsortedS[rL * blockSize+4]>alreadySorted[rR*blockSize+2]-alreadySorted[rR*blockSize+4]){
                                        comparisonInd=1;
                                    }
                                }
                            }
                        }
                    }
                    if(comparisonInd==1){
                        for(counter=0;counter<blockSize;counter++){
                            help[r*blockSize+counter]=unsortedS[rL*blockSize+counter];
                        }
                        rL++;
                    }
                    else{
                        for(counter=0;counter<blockSize;counter++){
                            help[r*blockSize+counter]=alreadySorted[rR*blockSize+counter];
                        }
                        rR++;
                    }
                }
            }
            r++;
        }
        for(r=0;r<length;r++){
            for(counter=0;counter<blockSize;counter++){
                unsortedS[r*blockSize+counter]=help[r*blockSize+counter];
            }
        }
        delete[] help;
    }
    return 0;
}


//*************************************************************************************
//***************************      CreateMemObjects
//*************************************************************************************
//
//
//  The function turns *myGraph into appropriate sequences that are kept in the global
//  memory of the graphic card. After the execution of this function the global memory
//  will contain the following sequences:
//
//            memObjects[0] - VERTICES
//  Each vertex is represented by 5 integers. They are name, label, status, firstEdge,
//  and relabelAttempt.
//  name is an integer from 0 to numberOfVertices-1.
//  label is the number that represents the quality of the water that is present at the
//  vertex.
//  status can be 1 (active) or 0 (inactive).
//  firstEdge is the position in the sequence of edges (i.e. memObjects[1]) of the first
//  edge originating from the vertex.
//  relabelAttempt is a storage location that will store the possible new label for
//  triggered vertex during one algorithm cycle.
//
//            memObjects[1] - EDGES
//  Each edge is represented by 8 integers. They are start, end, remainingTime,
//  f2 (weight/cost), f1 (original), source, status, opposite.
//  The edges will be listed twice. Namely the edge between vertices 0 and 1 and the edge
//  between 1 and 0 will be listed separately, and at all times all the parameters except
//  for start, end, and opposite will be the same. In the two listings the start and end
//  will be reversed and the variable opposite will contain the location of the opposite
//  edge to make it more convenient for program to make sure both sequences have all other
//  variables identical.
//  status of an edge can be 0 (passive); 1 (used); 2 (just used); 3 (active); 4 (active
//  and activated in this cycle).
//  The names of some variables suggest that sometimes edges with label 3 are referred to
//  as "ready" and edges with label 4 are referred to as "active."
//
//            memObjects[2] - NUMBER OF PARAMETERS (61, but some are unused)
//
//            memObjects[3] - SEQUENCE OF PARAMETERS
//  parameters[0]:  number of vertices
//  parameters[1]:  number of edges
//  parameters[2]:  number of integers for each vertex (5)
//  parameters[3]:  number of integers for each edge (8)
//  parameters[4]:  temporary storage of the length of a sequence during parallel
//                  cleaning.
//  parameters[6]:  output of terminalConditionCheck: -1 if still searching for path, or
//                  the last edge on the path if the path is found
//  parameters[7]:  dimension (or 1/2 of the bound of the degree of graph)
//  parameters[8]:  cost of the last edge on the shortest path (once the path is found)
//  parameters[9]:  number of active vertices
//  parameters[10]: number of active edges
//  parameters[11]: bound on the number of active vertices
//  parameters[12]: bound on the number of active edges
//  parameters[13]: bound on the number of elements of the auxiliary sequence
//  parameters[14]: has value 0 or 1 and is used by parallel merge sort to decide which
//                  sequence to sort. If the value is 0, then activeVertices are
//                  sorted; if the value is 1 then activeEdges are sorted.
//  parameters[17]: number of threads (usually called nThreads or numberOfThreads)
//  parameters[20]: removalSign (sign that means that the location is empty in the
//                  sequences of active edges, active vertices, and phantoms.
//                  The current value is -7.
//  parameters[27]: This parameter is useful in recovering the shortest path. It is set
//                  only when the terminalCondition is satisfied and  contains the label
//                  of the last vertex of the path (i.e. the label of the element in B
//                  in which the path terminates).
//  parameters[28]: number of elements in B
//  parameters[29]: number of phantoms
//  parameters[30]: bound on the number of phantoms
//  parameters[31]: number of integers for each phantom (5)
//  parameters[32]: the location in activeVertices after the last triggered vertex
//  parameters[33]: starting active edge for analysis in Step 1 of each algorithm cycle
//  parameters[34]: the location in activeEdges after the last triggered edge
//  parameters[35]: starting triggered vertex for analysis in Step 2of each algorithm
//                  cycle
//  parameters[36]: the location in the sequence of phantoms after the last phantom edge
//  parameters[37]: a number that denotes the empty slot in the storage for relabelAttempt
//  parameters[38]: used in parallel merge sort. Denotes the size of blocks that are
//                  already sorted
//  parameters[39]: used in parallel merge sort and in spotting the first empty slot
//                  after the cleaning
//                  When used in merge sort it has value of 0 or 1. The value denotes
//                  the destination of the merged sequence.
//                      - if it is 0 then the current sequence after the merge is
//                          located in auxiliary sequence
//                      - if it is 1 then the auxiliarry sequence is sorted, merged,
//                          and then transfered to the real sequence
//                  When used in cleaning stage it contains the location of the
//                  first empty slot in the sequence.
//  parameters[40]: the size of the memory of the sequence of active vertices (larger
//                  than the actual sequence could be)
//  parameters[41]: the size of the memory of the sequence of active edges (larger
//                  than the actual sequence could be)
//  parameters[42]: debugging mode
//  parameters[43]: if set to 1 then CPU is used instead of GPU
//
//            memObjects[4] - ACTIVE VERTICES
//  The sequence contains the locations of active vertices in the sequence of all vertices
//  (i.e. memObjects[0]). Their names can be obtained as the quotient
//      location / numberOfIntegersForEachVertex.
//
//            memObjects[5] - ACTIVE EDGES
//  The sequence contains the locations of active edges in the sequence of all edges
//  (i.e. memObjects[1]).
//  Their names can be obtained as location/numberOfIntegersForEachEdge.
//
//            memObjects[6] - AUXILIARY SEQUENCE
//  big sequence capabel of holding copies of each of the sequences of active vertices,
//  active edges, or phantoms. It is used during parallel mergeSort or parallel cleaning.
//
//
//            memObjects[7] - PHANTOMS
//  The sequence contains the phantoms. Each phantom is representing a phantom edge
//  and we use 5 terms in the sequence to store the following:
//  startPoint, endPoint, labelOfTheSource, remainingTime, weightOfTheEdge
//
//            memObjects[8] - B
//  The sequence of elements of B. Contains the location indices of elements of B. To
//  get the names of elements of B these locations should be divided by
//  numberOfIntegersForEachVertex.


myint CreateMemObjects(cl_context context, cl_mem *memObjects,
                       Graph *myGraph, myint *parameters, myint *numberOfParameters, myint numMemObjects, myint auxNumber)
{
    
    for(myint i=0;i<numMemObjects;i++){
        if(memObjects[i]!=0){
            clReleaseMemObject(memObjects[i]);
        }
    }
    
    VtxIt begV,endV;
    AdjIt begE,endE;
    EdgeDesc edge1;
    myint boundOnElementsInB=20;
    myint defaultNumberForEmpty=parameters[20];
        time_t tBeg,tEnd;
    tie(begV,endV)=vertices(*myGraph);
    myint numVertices=0;
    while(begV!=endV){
        begV++;
        numVertices++;
    }
    myint numEdges=0;
    EdgeIt edgeBeg,edgeEnd;
    for(tie(edgeBeg,edgeEnd)=edges(*myGraph);edgeBeg!=edgeEnd;++edgeBeg){
        numEdges++;
    }
    
    myint *verticesList;
    myint *edgesList;
    myint numberOfIntegersForEachVertex=parameters[2];
    myint numberOfIntegersForEachEdge=parameters[3];
    
    myint vCounterB, eCounterB;
    myint nThreadsB=parameters[17]+10;
    vCounterB=numVertices+nThreadsB;
    
    eCounterB=numEdges+nThreadsB;
    

    verticesList=new myint[vCounterB* numberOfIntegersForEachVertex];

    edgesList=new myint[2*eCounterB*numberOfIntegersForEachEdge];

    myint *elementsOfBList;
    elementsOfBList=new myint[boundOnElementsInB];
    

    myint numberOfElementsInB=0;
    parameters[0]=numVertices;
    parameters[1]=numEdges;
    parameters[40]=vCounterB;
    parameters[41]=eCounterB;
    tie(begV,endV)=vertices(*myGraph);
    myint vCounter=0, eCounter=0,eCounter2=0,eCounter3=0,totalECounter=0, totalVCounter=0;
    
    myint help1, help2;
    myint tempName,tempName2;
    myint tempBelongsToB;
    myint boundOnTheDegreeOfGraph=parameters[7] * 2;

    while(begV!=endV){
        tempName=get(&VertexProperties::vName,*myGraph,*begV);
        verticesList[vCounter]=tempName;
        vCounter++;
        verticesList[vCounter]=get(&VertexProperties::vCost,*myGraph,*begV);
        vCounter++;
        help1=(myint)(get(&VertexProperties::belongsToC,*myGraph,*begV));
        tempBelongsToB=(myint)(get(&VertexProperties::belongsToD,*myGraph,*begV));
        if(tempBelongsToB==2){
            tempBelongsToB=0;
            elementsOfBList[numberOfElementsInB]=tempName*numberOfIntegersForEachVertex;
            numberOfElementsInB++;
        }
        help1=2*help1+tempBelongsToB;
        help1=2*help1+(myint)(get(&VertexProperties::vState,*myGraph,*begV));
        verticesList[vCounter]=help1;
        vCounter++;
        verticesList[vCounter]=eCounter;
        vCounter++;
        
        verticesList[vCounter]=parameters[37];
        vCounter++;
        
        for(tie(begE,endE)=adjacent_vertices(*begV,*myGraph);begE!=endE;++begE){
            edge1=edge(*begV,*begE,*myGraph).first;
            tempName2=get(&VertexProperties::vName,*myGraph,*begE);
            eCounter3=eCounter;
            edgesList[eCounter]=tempName;
            eCounter++;
            edgesList[eCounter]=tempName2;
            eCounter++;
            edgesList[eCounter]=get(&EdgeProperties::eTime,*myGraph,edge1);
            eCounter++;
            edgesList[eCounter]=get(&EdgeProperties::eCost,*myGraph,edge1);
            eCounter++;
            edgesList[eCounter]=get(&EdgeProperties::eTime,*myGraph,edge1);
            eCounter++;
            edgesList[eCounter]=get(&EdgeProperties::eSource,*myGraph,edge1);
            eCounter++;
            edgesList[eCounter]=(myint)get(&EdgeProperties::eState,*myGraph,edge1);
            if(tempName<tempName2){
                put(&EdgeProperties::edgeId1,*myGraph,edge1,eCounter3);
            }
            else{
                put(&EdgeProperties::edgeId2,*myGraph,edge1,eCounter3);
            }
            
            
            eCounter++;
            eCounter++;
            
        }
        
        
        
        begV++;
    }

    
    totalECounter=eCounter;

    

    
    for(tie(edgeBeg,edgeEnd)=edges(*myGraph);edgeBeg!=edgeEnd;++edgeBeg){
        eCounter2=get(&EdgeProperties::edgeId1,*myGraph,*edgeBeg);
        eCounter3=get(&EdgeProperties::edgeId2,*myGraph,*edgeBeg);
        edgesList[eCounter2+numberOfIntegersForEachEdge-1]=eCounter3;
        edgesList[eCounter3+numberOfIntegersForEachEdge-1]=eCounter2;
        
    }
    
    totalVCounter=vCounter;
    
    
    myint totalVCounterB, totalECounterB;
    totalVCounterB=totalVCounter+ nThreadsB * numberOfIntegersForEachVertex;
    totalECounterB=totalECounter+ 2*nThreadsB *numberOfIntegersForEachEdge;
    
 
    memObjects[0] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * totalVCounterB , verticesList, NULL);

    memObjects[1] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * totalECounterB , edgesList, NULL);
    
    
    myint *activeVertices, *activeEdges, *auxSequence, *phantomSequence;
    myint auxNumHelpV=numVertices+auxNumber*boundOnTheDegreeOfGraph;
    myint auxNumHelpE=numEdges+auxNumber*boundOnTheDegreeOfGraph + 3* parameters[17];
    myint numberOfPhantoms=4 * parameters[31] * auxNumHelpE;
    
    if(numberOfPhantoms<10){
        numberOfPhantoms=10;
    }
    
    
    parameters[11]=auxNumHelpV;
    parameters[12]=auxNumHelpE;
    parameters[13]=numberOfPhantoms;
    memObjects[2] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) , numberOfParameters, NULL);
    
    
    activeVertices=new myint[auxNumHelpV];
    activeEdges=new myint[auxNumHelpE];
    phantomSequence=new myint[numberOfPhantoms];
    
    
    
    tie(begV,endV)=vertices(*myGraph);
    myint vCount1=0,vCount2=0,eCount1=0,eCount2=0;
    while(begV!=endV){
        if(get(&VertexProperties::vState,*myGraph,*begV)==1){
            activeVertices[vCount1]=vCount2;
            vCount1++;
        };
        
        
        
        for(tie(begE,endE)=adjacent_vertices(*begV,*myGraph);begE!=endE;++begE){
            edge1=edge(*begV,*begE,*myGraph).first;
            tempName2=get(&VertexProperties::vName,*myGraph,*begE);
            if( (get(&VertexProperties::vState,*myGraph,*begV)==1) || (get(&VertexProperties::vState,*myGraph,*begE)==1)){
                activeEdges[eCount1]=eCount2;
                eCount1++;
            }
            eCount2+=numberOfIntegersForEachEdge;
        }
        
        
        vCount2+=numberOfIntegersForEachVertex;
        begV++;
    }
    
    
    mergeSortDec(activeVertices,vCount1);
    mergeSortDec(activeEdges,eCount1);
    
    
    for(myint j=vCount1;j<auxNumHelpV;j++){
        activeVertices[j]=defaultNumberForEmpty;
    }
    
    parameters[9]=vCount1;
    parameters[10]=eCount1;
    parameters[28]=numberOfElementsInB;
    parameters[30]=numberOfPhantoms;
    
    
    
    for(myint j=0;j<numberOfPhantoms;j++){
        phantomSequence[j]=defaultNumberForEmpty;
    }
    
    
    parameters[29]=0;
    
    
    
    memObjects[3] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * (*numberOfParameters) , parameters, NULL);
   
    
    memObjects[4] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * auxNumHelpV , activeVertices, NULL);
 
    for(myint j=eCount1;j<auxNumHelpE;j++){
        activeEdges[j]=defaultNumberForEmpty;
    }
    memObjects[5] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * auxNumHelpE , activeEdges, NULL);
 

    for(myint j=0;j<auxNumHelpE;j++){
        activeEdges[j]=defaultNumberForEmpty;
    }
    
    
 

    
    memObjects[6] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * numberOfPhantoms , phantomSequence, NULL);

 

    
    
    memObjects[7] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * numberOfPhantoms , phantomSequence, NULL);
    
    memObjects[8] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                   sizeof(myint) * numberOfElementsInB , elementsOfBList, NULL);
    for(mysint j=0;j<numMemObjects;j++){
        if (memObjects[j] == NULL)
        {
            cerr << "Error creating memory objects.  " << j<<endl;
            return 0;
        }
    }
    
    
    delete[] edgesList;
    delete[] verticesList;
    delete[] activeVertices;
    delete[] activeEdges;
    delete[] elementsOfBList;
    delete[] phantomSequence;
    
    
    return auxNumHelpE;
}

void Cleanup(cl_context context, cl_command_queue commandQueue,
             cl_program program, cl_kernel *kernfls, myint numKernels, cl_mem *memObjects, myint numMemObjects)
{
    for (myint i = 0; i < numMemObjects; i++)
    {
        if (memObjects[i] != 0)
            clReleaseMemObject(memObjects[i]);
    }
    if (commandQueue != 0)
        clReleaseCommandQueue(commandQueue);
    
    for(myint k=0;k<numKernels;k++){
        if(kernfls[k]!=0){
            clReleaseKernel(kernfls[k]);
        }
    }

    
    
    
    if (program != 0)
        clReleaseProgram(program);
    
    if (context != 0)
        clReleaseContext(context);
    
}










myint treatError(cl_int errNum, cl_context context, cl_command_queue commandQueue,
                 cl_program program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint numMemObjects){
    std::cerr << "Error queuing kernel for execution. Kernel 1! " << " "<<errNum<< std::endl;
    if(errNum==CL_INVALID_PROGRAM_EXECUTABLE){cout<<"CL_INVALID_PROGRAM_EXECUTABLE"<<endl;}
    if(errNum==CL_INVALID_COMMAND_QUEUE){cout<<"CL_INVALID_COMMAND_QUEUE"<<endl;}
    if(errNum==CL_INVALID_KERNEL){cout<<"CL_INVALID_KERNEL"<<endl;}
    if(errNum==CL_INVALID_CONTEXT){cout<<"CL_INVALID_CONTEXT"<<endl;}
    if(errNum==CL_INVALID_KERNEL_ARGS){cout<<"CL_INVALID_KERNEL_ARGS"<<endl;}
    if(errNum==CL_INVALID_WORK_DIMENSION){cout<<"CL_INVALID_WORK_DIMENSION"<<endl;}
    if(errNum==CL_INVALID_GLOBAL_WORK_SIZE){cout<<"CL_INVALID_GLOBAL_WORK_SIZE"<<endl;}
    if(errNum==CL_INVALID_GLOBAL_OFFSET){cout<<"CL_INVALID_GLOBAL_OFFSET"<<endl;}
    if(errNum==CL_INVALID_WORK_GROUP_SIZE){cout<<"CL_INVALID_WORK_GROUP_SIZE"<<endl;}
    if(errNum==CL_INVALID_WORK_ITEM_SIZE){cout<<"CL_INVALID_WORK_ITEM_SIZE"<<endl;}
    if(errNum==CL_MISALIGNED_SUB_BUFFER_OFFSET){cout<<"CL_MISALIGNED_SUB_BUFFER_OFFSET"<<endl;}
    if(errNum==CL_OUT_OF_RESOURCES){cout<<"CL_OUT_OF_RESOURCES"<<endl;}
    if(errNum==CL_MEM_OBJECT_ALLOCATION_FAILURE){cout<<"CL_MEM_OBJECT_ALLOCATION_FAILURE"<<endl;}
    if(errNum==CL_INVALID_EVENT_WAIT_LIST){cout<<"CL_INVALID_EVENT_WAIT_LIST"<<endl;}
    if(errNum==CL_OUT_OF_HOST_MEMORY){cout<<"CL_OUT_OF_HOST_MEMORY"<<endl;}
    
    
    Cleanup(context, commandQueue, program, kernfls, numKernels, memObjects, numMemObjects);
    return 1;
}
mysint setMainKernelArguments(cl_kernel* kernel, cl_mem* memObjects, myint numMemObjects){
    cl_int errNum;
 
    for(mysint i=0;i<numMemObjects;i++){
        errNum = clSetKernelArg(*kernel, i, sizeof(cl_mem), &memObjects[i]);
        if (errNum != CL_SUCCESS)
        {
            cerr << "Error setting kernel argument " << i<<" "<<errNum<<endl;
            if(errNum==CL_INVALID_KERNEL){
                cout<<"CL_INVALID_KERNEL"<<endl;
            }
            if(errNum==CL_INVALID_ARG_INDEX){
                cout<<"CL_INVALID_ARG_INDEX"<<endl;
            }
            if(errNum==CL_INVALID_ARG_VALUE){
                cout<<"CL_INVALID_ARG_VALUE"<<endl;
            }
            if(errNum==CL_INVALID_MEM_OBJECT){
                cout<<"CL_INVALID_MEM_OBJECT"<<endl;
            }
            if(errNum==CL_INVALID_SAMPLER){
                cout<<"CL_INVALID_SAMPLER"<<endl;
            }
            if(errNum==CL_INVALID_ARG_SIZE){
                cout<<"CL_INVALID_ARG_SIZE"<<endl;
            }
            if(errNum==CL_OUT_OF_RESOURCES){
                cout<<"CL_OUT_OF_RESOURCES"<<endl;
            }
            if(errNum==CL_OUT_OF_HOST_MEMORY){
                cout<<"CL_OUT_OF_HOST_MEMORY"<<endl;
            }
        }
    }
     return 1;
}
//*************************************************************************************
//***************************      changeParameters
//*************************************************************************************
//
//
//  The function is used after the sequence of parameters is manipulated by CPU. The function updates memObjects[3] in the global memory.

mysint changeParameters(cl_context *context, 	cl_command_queue *commandqueue, cl_kernel *kernel,  cl_mem* memObjects, myint *parameters, myint numberOfParameters){

    cl_int errNum;
    cl_event eventE;
    errNum = clEnqueueWriteBuffer (*commandqueue, memObjects[3], CL_TRUE, 0, sizeof(myint) * numberOfParameters, parameters,0,NULL, &eventE);
    clWaitForEvents(1, &eventE);
    
    
    return 1;
}

//*************************************************************************************
//***************************      checkTerminalCondition
//*************************************************************************************
//
//
//  The function checks whether the algorithm has finished. The algorithm should stop
//  running once there is an element in B that is triggered
//  or activated, or if there are no more active edges.


myint checkTerminalCondition(cl_context *context, cl_command_queue *commandQueue,
                       cl_program *program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint* parameters, myint* numParameters, size_t sl, size_t pws, myint numMemObjects){
    myint terminalCondition=-1;
    cl_int errNum;
    cl_event eventE;
    size_t globalWorkSizeP[1]={ pws };
    size_t localWorkSizeP[1]={ pws };
    
 
    
    changeParameters(context,commandQueue, kernfls, memObjects, parameters, *numParameters);
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[0], 1, NULL,
                           globalWorkSizeP, localWorkSizeP,
                           0, NULL, &eventE);
     
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                 0,  (*numParameters) * sizeof(myint), parameters,
                                 0, NULL, &eventE);
    
    clWaitForEvents(1, &eventE);
    
    if(parameters[6]!=-1){
        terminalCondition=parameters[6];
    }
    if(terminalCondition==-1){
        if(parameters[10]==0){
            terminalCondition=-2;
        }
    }

    return terminalCondition;
}
//*************************************************************************************
//***************************      mergeSortDecParallel
//*************************************************************************************
//
//
//  The procedure sorts the sequence in decreasing order using parallel merge sort.
//  If parameters[14]=1, then the sequence that will be sorted is
//      the subsequence of activeVertices that starts at position parameters[32] and has length parameters[4]
//  If parameters[14]=2, then the sequence to be sorted is
//      the subsequence of activeEdges that starts at position parameters[34] and has length parameters[4]
myint mergeSortDecParallel(cl_context *context, cl_command_queue *commandQueue,
                           cl_program *program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint* parameters, myint* numParameters, size_t sl, size_t pws, myint numMemObjects){
    
    size_t globalWorkSizeP[1]={ pws };
    size_t localWorkSizeP[1]={ pws };
    size_t globalWorkSizePS[1]={ pws };
    size_t globalWorkSizePQ[1]={ pws };
    size_t localWorkSizePS[1]={ pws };
    
    cl_event eventE;
    cl_int errNum;
    myint sz1=0;
    
    myint blockSizePrev=1;
    myint blockSizeNew;
    myint cleaningLength=parameters[4];
    myint direction=0;
    
    myint kernCopyNum,  kernStep1MergeS, kernStep2MergeS, kernIncBSACD;
    myint  decisionNK;
    if(parameters[14]==1){
        kernCopyNum=10;
        kernStep1MergeS=12;
        kernStep2MergeS=13;
    }
    else{

        kernCopyNum=11;

        kernStep1MergeS=17;
        kernStep2MergeS=18;
    }
    
    parameters[38]=blockSizePrev;
    parameters[39]=0;
    
    kernIncBSACD=28;
    
    changeParameters(context, commandQueue, kernfls+kernCopyNum, memObjects, parameters, *numParameters);
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    // Each cycle of the following while loop assumes that
    // the sequences of length blockSizePrev are already sorted.
    // In each cycle two such consecutive sequences are merged.
    // The merging procedure consists of two stages.
    //      Stage 1.    For each element we calculate the location that it will have in the new sequence
    //                  The calculated locations are stored in remote locations of auxiliarySequence
    //      Stage 2.    Each element is copied to the correct location.
    //                  Unfortunately, the copying has to be performed to a new sequence, otherwise some locations may get erased
    //                  before their value is copied to the correct place.
    //                  The other sequence used is auxiliarySequence (but not remote locations of auxiliarySequence)
    //                  In order to improve performance we maintain another parameter called direction (parameter[39]) that is either 0 or 1
    //                  If it is 0, then the copying is done from the sequence to auxiliarySequence.
    //                  If it is 1, then the copying is done the other way round.
    while(blockSizePrev<cleaningLength){
        
        blockSizeNew= 2 * blockSizePrev;

        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep1MergeS], 1, NULL,
                                      globalWorkSizePQ, localWorkSizeP,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
        
        
        
        direction++;
        if(direction>1){
            direction=0;
        }
        
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep2MergeS], 1, NULL,
                                      globalWorkSizePQ, localWorkSizeP,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernIncBSACD], 1, NULL,
                                      globalWorkSizePS, localWorkSizePS,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
        
        blockSizePrev=blockSizeNew;
    }
    
    //  If the direction parameter is equal to 1 then the last copying happened from the original sequence to the auxiliarySequence
    //  and the auxiliarySequence has the final sorted sequence. It needs to be copied back to the original locations.
    if(direction==1){
        sz1= cleaningLength;
        if((sz1%pws!=0)||(sz1==0)){
            sz1=(sz1/pws+1)*pws;
        }
        globalWorkSizeP[0]=sz1;

        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernCopyNum], 1, NULL,
                                      globalWorkSizeP, localWorkSizeP,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
    }
    return 1;
    
}

//*************************************************************************************
//***************************      mergeSortDecBlockParallel
//*************************************************************************************
//
//
//  When soriting the sequence of phatnoms one should use mergeSortDecBlockParallel instead of mergeSortDecParallel.
//  The differences are analogous to the ones between mergeSortDecBlock and mergeSortDec.
myint mergeSortDecBlockParallel(cl_context *context, cl_command_queue *commandQueue,
                                cl_program *program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint* parameters, myint* numParameters, size_t sl, size_t pws, myint numMemObjects){
    
    size_t globalWorkSizeP[1]={ pws };
    size_t localWorkSizeP[1]={ pws };
    size_t globalWorkSizePS[1]={ pws };
    size_t globalWorkSizePQ[1]={ pws };
    size_t localWorkSizePS[1]={ pws };
    
    cl_event eventE;
    cl_int errNum;
    myint sz1=0;
    
    myint blockSizePrev=1;
    myint blockSizeNew;
    myint cleaningLength=parameters[4];
    myint direction=0;
    
    myint boundElements, kernCopyNum,  kernStep1MergeS, kernStep2MergeS, kernIncBSACD;
    myint  decisionNK;
    
    boundElements= parameters[30];
    kernCopyNum=24;
    kernStep1MergeS=22;
    kernStep2MergeS=23;

    
    parameters[38]=blockSizePrev;
    parameters[39]=0;
    
    kernIncBSACD=28;
    
    changeParameters(context, commandQueue, kernfls+kernStep1MergeS, memObjects, parameters, *numParameters);
    
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    while(blockSizePrev<cleaningLength){
        
        blockSizeNew= 2 * blockSizePrev;
        

        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep1MergeS], 1, NULL,
                                      globalWorkSizePQ, localWorkSizeP,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep2MergeS], 1, NULL,
                                      globalWorkSizePQ, localWorkSizeP,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
        
        direction++;
        if(direction>1){
            direction=0;
        }
        
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernIncBSACD], 1, NULL,
                                      globalWorkSizePS, localWorkSizePS,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        blockSizePrev=blockSizeNew;
    }
    if(direction==1){
        sz1= cleaningLength;
        if((sz1%pws!=0)||(sz1==0)){
            sz1=(sz1/pws+1)*pws;
        }
        globalWorkSizeP[0]=sz1;
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernCopyNum], 1, NULL,
                                      globalWorkSizeP, localWorkSizeP,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
    }
    return 1;
    
}

//*************************************************************************************
//***************************      cleanSequence_parallel
//*************************************************************************************
//
//
//  "Cleaning" the sequence means sorting it in decreasing order, removing any repetitions of terms, and removing any terms that
//  are equal to removalSign (i.e. -7)
//  If parameters[14]=1 the sequence that is cleaned is the subsequence of activeVertices from position parameters[32] and of length parameters[4]
//  If parameters[14]=2 the sequence that is cleaned is the subsequence of activeEdges from position parameters[34] and of length parameters[4]

myint cleanSequence_parallel(cl_context *context, cl_command_queue *commandQueue,
                                   cl_program *program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint* parameters, myint* numParameters, size_t sl, size_t pws, myint numMemObjects){
    size_t localWorkSizeP[1]={ pws };
    size_t globalWorkSizePS[1]={ pws };
    size_t globalWorkSizePQ[1]={ pws };
    size_t localWorkSizePS[1]={ pws };
    
    cl_event eventE;
    cl_int errNum;
    myint sz1=0;
    
    myint blockSizePrev=1;
    myint blockSizeNew;

    myint direction=0;
    
    myint  kernStep1Clean, kernStep2Clean, kernSpotRemSign, decisionNK;
    if(parameters[14]==1){
        kernStep1Clean=14;
        kernStep2Clean=15;
        kernSpotRemSign=16;
        parameters[4]+=parameters[32]-parameters[9];
        parameters[32]=parameters[9];
    }
     else{
         kernStep1Clean=19;
         kernStep2Clean=20;
         kernSpotRemSign=21;
         parameters[4]+=parameters[34]-parameters[10];
         parameters[34]=parameters[10];

     }
    myint cleaningLength=parameters[4];

    
    parameters[38]=blockSizePrev;
    parameters[39]=0;
    

    
    mergeSortDecParallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
    
    
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernSpotRemSign], 1, NULL,
                                  globalWorkSizePS, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                 0,  (*numParameters) * sizeof(myint), parameters,
                                 0, NULL, &eventE);
    
    clWaitForEvents(1, &eventE);
    
    if(parameters[14]==1){
        if (parameters[39]<parameters[32]+cleaningLength){
            cleaningLength=parameters[39]-parameters[32];
            parameters[4]=cleaningLength;
        }
    }
    else{
        if (parameters[39]<parameters[34]+cleaningLength){
            cleaningLength=parameters[39]-parameters[34];
            parameters[4]=cleaningLength;
        }
    }
    
    changeParameters(context, commandQueue, kernfls+kernStep1Clean, memObjects, parameters, *numParameters);
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep1Clean], 1, NULL,
                                  globalWorkSizePQ, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep2Clean], 1, NULL,
                                  globalWorkSizePQ, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    
    mergeSortDecParallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
    
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernSpotRemSign], 1, NULL,
                                  globalWorkSizePS, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                 0,  (*numParameters) * sizeof(myint), parameters,
                                 0, NULL, &eventE);
    
    clWaitForEvents(1, &eventE);
    if(parameters[14]==1){
        parameters[32]=parameters[39];
    }
    else{
        parameters[34]=parameters[39];
    }
 
    
    
    
    return 1;
    
}


//*************************************************************************************
//***************************      cleanPhantomSequence_parallel
//*************************************************************************************
//
//
//  Analogous to cleanSequence_parallel except that this one is applied to the sequence of phantoms which consists of blocks of length 5
//  and that are compared in a specific way described in the comment near the function mergeSortDecBlock.


myint cleanPhantomSequence_parallel(cl_context *context, cl_command_queue *commandQueue,
                                     cl_program *program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint* parameters, myint* numParameters, size_t sl, size_t pws, myint numMemObjects){
    size_t localWorkSizeP[1]={ pws };
    size_t globalWorkSizePS[1]={ pws };
    size_t globalWorkSizePQ[1]={ pws };
    size_t localWorkSizePS[1]={ pws };
    
    cl_event eventE;
    cl_int errNum;
    myint sz1=0;
    
    myint blockSizePrev=1;
    myint blockSizeNew;
    
    myint direction=0;
    
    
    myint  kernStep1Clean, kernStep2Clean, kernSpotRemSign, decisionNK;
    
    
    
    
    
    
    
    kernStep1Clean=25;
    kernStep2Clean=26;
    kernSpotRemSign=27;
    parameters[4]+=parameters[36]-parameters[29];
    parameters[36]=parameters[29];
    
    
    
    if(parameters[4]==0){return 1;}
    
    myint cleaningLength=parameters[4];
    
    
    parameters[38]=blockSizePrev;
    parameters[39]=0;
    
    
    

    
    mergeSortDecBlockParallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
    
    
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernSpotRemSign], 1, NULL,
                                  globalWorkSizePS, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                 0,  (*numParameters) * sizeof(myint), parameters,
                                 0, NULL, &eventE);
    
    clWaitForEvents(1, &eventE);
    
    
    if (parameters[39]<parameters[36]+cleaningLength){
        cleaningLength=parameters[39]-parameters[36];
        parameters[4]=cleaningLength;
    }
    
    changeParameters(context, commandQueue, kernfls+kernStep1Clean, memObjects, parameters, *numParameters);
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep1Clean], 1, NULL,
                                  globalWorkSizePQ, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernStep2Clean], 1, NULL,
                                  globalWorkSizePQ, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    
    mergeSortDecBlockParallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
    
    decisionNK=parameters[4];
    if((decisionNK==0)||(decisionNK%pws!=0)){
        decisionNK=( parameters[4]/ pws + 1)*pws;
    }
    globalWorkSizePQ[0]=decisionNK;
    
    errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[kernSpotRemSign], 1, NULL,
                                  globalWorkSizePS, localWorkSizePS,
                                  0, NULL, &eventE);
    if (errNum != CL_SUCCESS){
        treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    clWaitForEvents(1, &eventE);
    errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                 0,  (*numParameters) * sizeof(myint), parameters,
                                 0, NULL, &eventE);
    
    clWaitForEvents(1, &eventE);
    
    parameters[36]=parameters[39];
    
    
    
    
    return 1;
    
}





//*************************************************************************************
//***************************      cleanSequence
//*************************************************************************************
//
//
//  When sequences are relatively short, we use this function instead of the parallel one.
//  GPUs are not as fast as the theory suggest...
//
// This procedure assumes that from *storageStart to *cleanFrom the sequence is sorted in decreasing order
// and that it does not contain symbols equal to removalSign.
// The procedure will find all the terms from *cleanFrom to (*cleanFrom+cleaningLength) that are different than
// removal sign, sort them and condense. *cleanFrom will contain the new length of the sequence.
myint cleanSequence(myint *sequence, myint *sequenceLength, myint *storageStart, myint *cleanFrom, myint cleaningLength, myint removalSign){
    myint storageEnd=*cleanFrom;
    myint readPos=*cleanFrom;
    myint writePos;
    myint currentTerm;
    myint readEnd=*cleanFrom+cleaningLength;
    myint counter=0; myint counter2;
    myint ind1;
    
    myint *seq, *auxSeq;
    seq=sequence+*cleanFrom;
    

    mergeSortDec(seq,cleaningLength);
    
    
    auxSeq=new myint[readEnd- *storageStart];
    myint position1, position2, currentPosition, biggerTerm,increaseHappened;
    position1=*storageStart;
    position2=0; currentPosition=0;
    while((position1<*cleanFrom)||(position2<cleaningLength)){
        biggerTerm=-1;
        
        increaseHappened=0;
        if(position1<*cleanFrom){
            biggerTerm=sequence[position1];
        }
        if(seq[position2]==removalSign){
            position2=cleaningLength;
        }
        if(position2<cleaningLength){
            if(seq[position2]>biggerTerm){
                biggerTerm=seq[position2];
                position2++;
                increaseHappened=1;
            }
        }
        if(increaseHappened==0){
            position1++;
        }
        if(((currentPosition==0)||(auxSeq[currentPosition-1]>biggerTerm))&&(biggerTerm>0)){
            auxSeq[currentPosition]=biggerTerm;
            currentPosition++;
        }
        
        
    }
    
    *cleanFrom=currentPosition+(*storageStart);
    writePos=*storageStart;
    for(myint i=0;i<currentPosition;i++){
        sequence[writePos]=auxSeq[i];
        writePos++;
    }
    while(writePos<readEnd){
        sequence[writePos]=removalSign;
        writePos++;
    }
    
    delete [] auxSeq;
    
    return 1;
    
}

//*************************************************************************************
//***************************      cleanPhantomSequence
//*************************************************************************************
//
//
//  Analogous to cleanSequence except that this one is used for the sequence of phantoms.
myint cleanPhantomSequence(myint *sequence, myint *sequenceLength, myint *storageStart, myint *cleanFrom, myint cleaningLength, myint removalSign,myint blockSize ){
    
    // This procedure assumes that from *storageStart to *cleanFrom the sequence is sorted in decreasing order
    // and that it does not contain symbols equal to removalSign.
    // The procedure will find all the terms from *cleanFrom to (*cleanFrom+cleaningLength) that are different than
    // removal sign and put them in storage. The *cleanFrom will contain the value of the end of the updated
    // storage.
    
    
    myint storageEnd=*cleanFrom;
    myint readPos=*cleanFrom;
    myint writePos;
    myint currentTerm;
    myint readEnd=*cleanFrom+cleaningLength;
    myint counter=0; myint counter2;
    myint ind1;
    myint cnter;
    myint *seq, *auxSeq;
    myint comparisonInd;
    
    
    seq=sequence+(*cleanFrom)*blockSize;
    
    
    
    
    mergeSortDecBlock(seq,cleaningLength,blockSize);
    
    
    auxSeq=new myint[(readEnd- *storageStart)*blockSize];
    myint position1, position2, currentPosition, increaseHappened;
    myint *biggerTerm;
    biggerTerm=new myint[blockSize];
    position1=*storageStart;
    position2=0; currentPosition=0;
    while((position1<*cleanFrom)||(position2<cleaningLength)){
        biggerTerm[1]=-1;
        
        increaseHappened=0;
        if(position1<*cleanFrom){
            for(cnter=0;cnter<blockSize;cnter++){
                biggerTerm[cnter]=sequence[position1*blockSize+cnter];
            }
        }
        if(seq[position2*blockSize+1]==removalSign){
            position2=cleaningLength;
        }
        if(position2<cleaningLength){
            comparisonInd=0;
            //  comparisonInd=seq[position2]>biggerTerm then comparisonInd=1
            if(seq[position2*blockSize+1]>biggerTerm[1]){
                comparisonInd=1;
            }
            else{
                if(seq[position2*blockSize+1]==biggerTerm[1]){
                    if(seq[position2*blockSize+3]>biggerTerm[3]){
                        comparisonInd=1;
                    }
                    else{
                        if(seq[position2*blockSize+3]==biggerTerm[3]){
                            if(seq[position2*blockSize+2]-seq[position2*blockSize+4]<biggerTerm[2]-biggerTerm[4]){
                                for(cnter=0;cnter<blockSize;cnter++){
                                    seq[position2*blockSize+cnter]=biggerTerm[cnter];
                                }
                            }
                            comparisonInd=1;
                        }
                    }
                }
            }
            if(comparisonInd==1){
                for(cnter=0;cnter<blockSize;cnter++){
                    biggerTerm[cnter]=seq[position2*blockSize+cnter];
                }
                position2++;
                increaseHappened=1;
            }
        }
        if(increaseHappened==0){
            position1++;
        }
        comparisonInd=0;
        if(currentPosition>0){
            //  comparisonInd=auxSeq[currentPosition-1]>biggerTerm then comparisonInd=1
            if(auxSeq[(currentPosition-1)*blockSize+1]>biggerTerm[1]){
                comparisonInd=1;
            }
            else{
                if(auxSeq[(currentPosition-1)*blockSize+1]==biggerTerm[1]){
                    if(auxSeq[(currentPosition-1)*blockSize+3]>biggerTerm[3]){
                        comparisonInd=1;
                    }
                    else{
                        if(auxSeq[(currentPosition-1)*blockSize+3]==biggerTerm[3]){
                            if(auxSeq[(currentPosition-1)*blockSize+2]-auxSeq[(currentPosition-1)*blockSize+4]<biggerTerm[2]-biggerTerm[4]){
                                for(cnter=0;cnter<blockSize;cnter++){
                                    auxSeq[(currentPosition-1)*blockSize+cnter]=biggerTerm[cnter];
                                }
                            }
                        }
                    }
                }
            }
            
        }
        if(((currentPosition==0)||(comparisonInd==1))&&(biggerTerm[1]>0)){
            for(cnter=0;cnter<blockSize;cnter++){
                auxSeq[currentPosition*blockSize+cnter]=biggerTerm[cnter];
            }
            currentPosition++;
        }
        
        
    }
    
    *cleanFrom=currentPosition+(*storageStart);
    writePos=*storageStart;
    for(myint i=0;i<currentPosition;i++){
        for(cnter=0;cnter<blockSize;cnter++){
            sequence[writePos*blockSize+cnter]=auxSeq[i*blockSize+cnter];
        }
        writePos++;
    }
    while(writePos<readEnd){
        for(cnter=0;cnter<blockSize;cnter++){
            sequence[writePos*blockSize+cnter]=removalSign;
        }
        writePos++;
    }
    
    delete [] auxSeq;
    delete[] biggerTerm;
    return 1;
    
}

//*************************************************************************************
//***************************      condenseSequence
//*************************************************************************************
//
//
//  A faster procedure then cleanSequence that can be used if the sequence is already sorted.
//


myint condenseSequence(myint *sequence, myint *phLength, myint removalSign){
    
    
    
    myint cnter;
    myint readLoc, writeLoc;
    readLoc=0;writeLoc=0;
    while(readLoc<*phLength){
        if(sequence[readLoc]!=removalSign){
            if(readLoc!=writeLoc){
                sequence[writeLoc]=sequence[readLoc];
            }
            writeLoc++;
            
        }
        readLoc++;
    }
    *phLength=writeLoc;
    while(writeLoc<readLoc){
        sequence[writeLoc]=removalSign;
        writeLoc++;
    }
    
    return 1;
    
}



//*************************************************************************************
//***************************      algCycle
//*************************************************************************************
//
//
//  Each cycle in the algorithm corresponds to one unit of time. The cycle consists of 9 steps.
//

myint algCycle(cl_context *context, cl_command_queue *commandQueue,
              cl_program *program, cl_kernel *kernfls, myint numKernels,cl_mem *memObjects, myint* parameters, myint* numParameters, size_t sl, size_t pws, myint numMemObjects, myint *activeVerticesSeq, myint *activeEdgesSeq, myint *phantomSeq){
    
    

    myint thresholdForParallel=40000;
    myint thresholdForParallel2=8000;
    size_t locSize1;
    size_t glSize1, glSizeTemp, glSizeTemp2, glSizeTemp3;
    cl_int errNum;
    cl_event eventE;
    size_t pws1;
    myint nThreads=parameters[17];
    myint removalSign=parameters[20];

    myint dimension=parameters[7];
    myint fCounter=0;
    myint nThreadsTemp, nThreadsTemp2, nThreadsTemp3;

    locSize1=pws;
    if((nThreads%pws!=0)||(nThreads==0)){
        nThreads=(nThreads/pws+1)*pws;
    }
    
    
    parameters[17]=nThreads;
    glSize1=nThreads;
    glSizeTemp=glSize1;
    size_t globalWorkSize[1] = { glSize1 };
    size_t localWorkSize[1] = { pws };
    size_t globalWorkSizeP[1]={ 1 };
    size_t localWorkSizeP[1]={ 1 };
    size_t globalWorkSizeTemp[1] = { glSizeTemp };
    
    

    parameters[24]=0;
 
    myint checkEnd;
    myint keeperOfNumberOfPhantoms, numberOfReadyEdges, numberOfPhantoms;

    myint paramSaver1, paramSaver2;
    changeParameters(context, commandQueue, kernfls+1, memObjects, parameters, *numParameters);
    
    errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                 0,  (*numParameters) * sizeof(myint), parameters,
                                 0, NULL, &eventE);
    
    clWaitForEvents(1, &eventE);
    
    //*********************************************************
    //*************     Step 1: Vertex triggering
    //*********************************************************
    //
    //  In this step we go over all active edges and decrease their
    //  time parameters by $1$. %Such edges are kept in the sequence
    //  of active edges. Once their parameter is decreased we change
    //  their status from ready to active.
    //  If for any edge the time parameter becomes $0$, the edge
    //  becomes {\em just used} and its destination triggered.
    //  We look at the destination of this edge and need to decide whether
    //  this destination will be triggered.
    //  When making such decision we can avoid triggering a vertex twice due
    //  to the fact that each edge appears twice in the sequence of active edges
    //  (the second appearance is when the endpoints are switched). Thus only
    //  if the destination of the water flow is actually the start point of the
    //  edge in consideration, we will perform an analysis to see whether it is
    //  going to be triggered.
    //
    //  If $Q$ is the destination and if the new water particle would increase the
    //  label of $Q$, then $Q$ is considered triggered and is added to the end of
    //  the sequence of active vertices. We are not modifying any vertices now
    //  because of the danger that several processing elements may be analyzing
    //  the same vertex.
    //
    //  To avoid the danger of two processing elements writing in the same location
    //  of the sequence of active vertices, we have to make sure that each processing
    //  element that runs concurrently has pre-specified location to write. This is
    //  accomplished by first specifying the number of threads in the separate variable
    //  nThreads}. Whenever kernels are executed in parallel we are using only nThreads
    //  processing elements. Each processing element has its id number which is used to
    //  determine the memory location to which it is allowed to write.
    //  The sequence of triggered vertices has to be cleaned after each parallel execution
    //  and at that point we take an additional step to ensure we don't list any of the
    //  vertices as triggered twice.


    
    numberOfReadyEdges=parameters[10];
    parameters[32]=parameters[9];
    nThreadsTemp=nThreads;
    glSizeTemp=glSize1;
    if(numberOfReadyEdges<nThreads){
        nThreadsTemp= (numberOfReadyEdges/ pws + 1) * pws;
        glSizeTemp=nThreadsTemp;
    }
    globalWorkSizeTemp[0]=glSizeTemp;
    glSizeTemp3=glSizeTemp;
    nThreadsTemp3=nThreadsTemp;
    while(fCounter<numberOfReadyEdges){
        parameters[33]=fCounter;
        changeParameters(context, commandQueue,kernfls+1, memObjects, parameters, *numParameters);
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[1], 1, NULL,
                                      globalWorkSizeTemp, localWorkSize,
                                      0, NULL, &eventE);

        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);

        errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                     0,  (*numParameters) * sizeof(myint), parameters,
                                     0, NULL, &eventE);
        
        clWaitForEvents(1, &eventE);
        if(nThreadsTemp>thresholdForParallel){
            parameters[14]=1;
            parameters[4]=nThreadsTemp;
            cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
        }
        else{
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[4], CL_TRUE, parameters[9]* sizeof(myint),  (parameters[32] + nThreadsTemp- parameters[9]) * sizeof(myint), activeVerticesSeq+parameters[9], 0, NULL, &eventE);
            
            clWaitForEvents(1, &eventE);
            paramSaver2=parameters[32];
            
            
            cleanSequence(activeVerticesSeq, parameters+11, parameters+9, parameters+32,nThreadsTemp,parameters[20] );
            
            clEnqueueWriteBuffer(*commandQueue, memObjects[4], CL_TRUE, parameters[9] * sizeof(myint), sizeof(myint) * (paramSaver2 + nThreadsTemp-parameters[9]), activeVerticesSeq + parameters[9],0,NULL, &eventE);
            clWaitForEvents(1, &eventE);
        }
        fCounter+=nThreadsTemp;
    }
    
    //*********************************************************
    //*************   Step 2: Analyzing triggered vertices
    //*********************************************************
    //
    //  For each triggered vertex $Q$ we look at all of its edges
    //  that are just used. We identify the largest possible label
    //  that can result from one of just used edges that starts
    //  from $Q$. That label will be stored in the sequence of
    //  vertices at the position reserved for temporary replacement
    //  label. The vertex is labeled as just triggered. If the vertex
    //  $Q$ is not active, this label will replace the current label
    //  of the vertex in one of the later steps. If the vertex $Q$
    //  is active, then this temporary label will be used later to
    //  construct an appropriate phantom edge.
    //  We are sure that different processing elements are not accessing
    //  the same vertex at the same time, because before this step we
    //  achieved the state in which there are no repetitions in the
    //  sequence of triggered vertices.

    
    fCounter=0;
    
    myint numberOfTriggeredVertices=parameters[32]-parameters[9];
    
    nThreadsTemp2=nThreadsTemp;
    glSizeTemp2=glSizeTemp;
    
    if(numberOfTriggeredVertices<nThreadsTemp){
        nThreadsTemp2= (numberOfTriggeredVertices/ pws + 1) * pws;
        glSizeTemp2=nThreadsTemp2;
    }
    
    globalWorkSizeTemp[0]=glSizeTemp2;
    

    while(fCounter<numberOfTriggeredVertices){
        parameters[35]=fCounter;
        
        
        changeParameters(context, commandQueue,kernfls+2, memObjects, parameters, *numParameters);
 
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[2], 1, NULL,
                                      globalWorkSizeTemp, localWorkSize,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
 
        fCounter+=nThreadsTemp2;
        
        
        
    }
    
    

    
    
    errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                 0,  (*numParameters) * sizeof(myint), parameters,
                                 0, NULL, &eventE);
    
    clWaitForEvents(1, &eventE);
    

    if(nThreadsTemp> thresholdForParallel){
        parameters[14]=1;
        parameters[4]=nThreadsTemp;
        cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
    }
    else{
        errNum = clEnqueueReadBuffer(*commandQueue, memObjects[4], CL_TRUE, parameters[9] * sizeof(myint),  (parameters[32] + nThreadsTemp - parameters[9]) * sizeof(myint), activeVerticesSeq + parameters[9], 0, NULL, &eventE);
        
        clWaitForEvents(1, &eventE);
        paramSaver2=parameters[32];
        
        
        cleanSequence(activeVerticesSeq, parameters+11, parameters+9, parameters+32,nThreadsTemp,parameters[20] );
        
        clEnqueueWriteBuffer(*commandQueue, memObjects[4], CL_TRUE, parameters[9] * sizeof(myint), sizeof(myint) * (paramSaver2 + nThreadsTemp-parameters[9]), activeVerticesSeq + parameters[9],0,NULL, &eventE);
        clWaitForEvents(1, &eventE);
    }
    
    
    //*********************************************************
    //*************   Step 3: Taking input from phantoms
    //*********************************************************
    //
    // The need to have this step separated from the previous ones
    //  is the current architecture of graphic cards that creates
    //  difficulties with dynamic memory locations. It is more
    //  efficient to keep phantom edges separate from the regular
    //  edges. The task is to look for all phantom edges and decrease
    //  their time parameters. If a phantom edge gets its time parameter
    //  equal to $0$, its destination is studied to see whether it
    //  should be added to the sequence of triggered edges. We calculate
    //  the new label that the vertex would receive through this phantom.
    //  We check whether this new label is higher than the currently known
    //  label and the temporary label from possibly previous triggering of
    //  the vertex.  The phantoms will not result in the concurrent writing
    //  to memory locations because each possible destination of a phantom
    //  could have only one edge that has time component equal to $0$.
    

    
    numberOfPhantoms=parameters[29];
    
    fCounter=0;

    
    
    nThreadsTemp=nThreads;
    glSizeTemp=glSize1;
    
    if(numberOfPhantoms<nThreads){
        nThreadsTemp= (numberOfPhantoms/ pws + 1) * pws;
        glSizeTemp=nThreadsTemp;
    }
    
    globalWorkSizeTemp[0]=glSizeTemp;
    

    
    while(fCounter<numberOfPhantoms){
        parameters[33]=fCounter;

        changeParameters(context, commandQueue, kernfls+3, memObjects, parameters, *numParameters);
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[3], 1, NULL,
                                      globalWorkSizeTemp, localWorkSize,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                     0,  (*numParameters) * sizeof(myint), parameters,
                                     0, NULL, &eventE);
        
        clWaitForEvents(1, &eventE);
        


        if(nThreadsTemp>thresholdForParallel){
            parameters[14]=1;
            parameters[4]=nThreadsTemp;
            cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
        }
        else{
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[4], CL_TRUE, parameters[9] * sizeof(myint),  (parameters[32] + nThreadsTemp-parameters[9]) * sizeof(myint), activeVerticesSeq + parameters[9], 0, NULL, &eventE);
            
            clWaitForEvents(1, &eventE);
            paramSaver2=parameters[32];
            
            
            cleanSequence(activeVerticesSeq, parameters+11, parameters+9, parameters+32,nThreadsTemp,parameters[20] );
            
            clEnqueueWriteBuffer(*commandQueue, memObjects[4], CL_TRUE, parameters[9] * sizeof(myint), sizeof(myint) * (paramSaver2 + nThreadsTemp-parameters[9]), activeVerticesSeq + parameters[9],0,NULL, &eventE);
            clWaitForEvents(1, &eventE);
            
        }

        
        fCounter+=nThreadsTemp;
    }
    
 
    //*********************************************************
    //*************   Step 4: Triggering edges
    //*********************************************************
    //
    //  Triggered vertices are analyzed using separate processing
    //  elements. A processing element analyzes the vertex $Q$
    //  in the following way.
    //  Each edge $j$ of $Q$ %that is not just used will be consider
    //  triggered if it can cause the other endpoint to get better
    //  label in future through $Q$. The edge $j$ is placed in
    //  the end of the sequence of active edges.
    


    fCounter=0;
    
    parameters[34]=parameters[10];
    numberOfTriggeredVertices=parameters[32]-parameters[9];
    
    nThreadsTemp=nThreadsTemp3;
    glSizeTemp=glSizeTemp3;
    
    if(numberOfTriggeredVertices<nThreadsTemp){
        nThreadsTemp= (numberOfTriggeredVertices/ pws + 1) * pws;
        glSizeTemp=nThreadsTemp;
    }
    
    globalWorkSizeTemp[0]=glSizeTemp;
    
    while(fCounter<numberOfTriggeredVertices){
        parameters[35]=fCounter;

        changeParameters(context,commandQueue, kernfls+4, memObjects, parameters, *numParameters);
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[4], 1, NULL,
                                      globalWorkSizeTemp, localWorkSize,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                     0,  (*numParameters) * sizeof(myint), parameters,
                                     0, NULL, &eventE);
        clWaitForEvents(1, &eventE);
        if(nThreadsTemp * dimension * 2>thresholdForParallel2){
            parameters[14]=2;
            parameters[4]=nThreadsTemp * dimension * 2;
            cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
        }
        else{
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[5], CL_TRUE, parameters[10] * sizeof(myint),  (parameters[34] + nThreadsTemp * dimension * 2 - parameters[10]) * sizeof(myint), activeEdgesSeq + parameters[10], 0, NULL, &eventE);
            
            clWaitForEvents(1, &eventE);
            paramSaver2=parameters[34];
            
            
            cleanSequence(activeEdgesSeq, parameters+12, parameters+10, parameters+34,nThreadsTemp * dimension * 2,parameters[20] );
            
            clEnqueueWriteBuffer(*commandQueue, memObjects[5], CL_TRUE, parameters[10] * sizeof(myint), sizeof(myint) * (paramSaver2 + nThreadsTemp*dimension*2 - parameters[10]), activeEdgesSeq + parameters[10],0,NULL, &eventE);
            clWaitForEvents(1, &eventE);

        }
        fCounter+=nThreadsTemp;
    }
 
    //*********************************************************
    //*************   Step 5: Treatment of triggered edges
    //*********************************************************
    //
    //  Each triggered edge will be analyzed with a dedicated
    //  processing element. %This analysis is performed in parallel.
    //  Consider a triggered edge $j$. We first identify its two
    //  endpoints. For the purposes of this step we will identify
    //  the endpoint with the larger label, call it the source, and
    //  denote by $S$. The other will be called the destination and
    //  denoted by $D$. In the end of the cycle, this vertex $S$ will
    //  become the source of the flow through $j$.
    //
    //  Notice that at least one of the endpoints is triggered.
    //  If only one endpoint is triggered, than we are sure that
    //  this triggered endpoint is the one that we designated as
    //  the source $S$.
    //
    //  We then look whether the source $S$ was active or inactive
    //  before it was triggered.
    //      Case 1: The source $S$ was inactive before triggering:
    //          There are several cases based on the prior status of $j$.
    //              If $j$ was passive, then it should become active
    //                  and no further analysis is necessary.
    //              If it was used or just used, then it should become
    //                  active and the time component should be restored
    //                  to the original one.
    //              Assume now that the edge $j$ was active. Based on
    //                  the knowledge that $S$ was inactive vertex we can
    //                  conclude that the source of $j$ was $D$. However
    //                  we know that the source of $j$ should be $S$ and
    //                  hence the time component of $j$ should be restored
    //                  to the backup value.
    //          Consequently, in the case that $S$ was inactive, regardless
    //          of what the status of $j$ was, we are sure its new status
    //          must be active and its time component can be restored to
    //          the original value. This restoration is not necessary in
    //          the case that $j$ was passive, however there is no harm in
    //          doing it.
    //          If the edge $j$ was not active before, then the edge $j$
    //          should be added to the list of active edges. This is
    //          accomplished by adding  only the opposite edge to the end
    //          of the list of triggered edges.
    //          If the edge $j$ was active before, then it should be removed
    //          from the list of triggered edges because all triggered edges
    //          will be merged into active edges. The edge $j$ already appears
    //          in the list of active edges and need not be added again.
    //      Case 2: The source $S$ was active before triggering:
    //          In this case we create phantom edges. Each such triggered edge
    //          generates four entries in the phantom sequence. The first one
    //          is the source, the second is the destination, the third is the
    //          label of the source (or the label stored in the temporary label
    //          slot, if higher), and the fourth is the original passage time
    //          through the edge $j$.
    //          The procedure described is realized through a loop. In each
    //          iteration, each thread may write in the sequence of phantoms.
    //          Thus, after each iteration the sequence of phantoms has to be
    //          cleaned to maintain its length.

    

    fCounter=0;
    
    
    myint numberOfTriggeredEdges=parameters[34]-parameters[10];
    numberOfPhantoms=parameters[29];
    parameters[36]=parameters[29];
    nThreadsTemp=nThreads;
    glSizeTemp=glSize1;
    if(numberOfTriggeredEdges<nThreads){
        nThreadsTemp= (numberOfTriggeredEdges/ pws + 1) * pws;
        glSizeTemp=nThreadsTemp;
    }
    globalWorkSizeTemp[0]=glSizeTemp;
    while(fCounter<numberOfTriggeredEdges){
        parameters[35]=fCounter;

        changeParameters(context, commandQueue,kernfls+5, memObjects, parameters, *numParameters);
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[5], 1, NULL,
                                      globalWorkSizeTemp, localWorkSize,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                     0,  (*numParameters) * sizeof(myint), parameters,
                                     0, NULL, &eventE);
        clWaitForEvents(1, &eventE);
        paramSaver2=(parameters[36]+ nThreadsTemp) * parameters[31];
        if(paramSaver2>thresholdForParallel2){
            parameters[4]=nThreadsTemp;
            cleanPhantomSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
        }
        else{
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[7], CL_TRUE,
                                         0,  paramSaver2 * sizeof(myint), phantomSeq,
                                         0, NULL, &eventE);
            
            clWaitForEvents(1, &eventE);
            keeperOfNumberOfPhantoms=parameters[29];
            parameters[29]=0;
            cleanPhantomSequence(phantomSeq, parameters+30, parameters+29,parameters+36,nThreadsTemp,parameters[20],parameters[31]);
            parameters[29]=parameters[36];
            clEnqueueWriteBuffer(*commandQueue, memObjects[7], CL_TRUE, 0, sizeof(myint) * paramSaver2, phantomSeq,0,NULL, &eventE);
            clWaitForEvents(1, &eventE);
        }
        fCounter+=nThreadsTemp;
    }
    
    changeParameters(context,commandQueue, kernfls+5, memObjects, parameters, *numParameters);
    
     paramSaver1=2*parameters[34];
    if(paramSaver1==0){
        paramSaver1=20;
    }
    if(paramSaver1>parameters[12]){
        paramSaver1=parameters[12];
    }

   
    
    myint cleaningLength1=(parameters[34]-parameters[10])*2;
    parameters[34]=parameters[10];
    if(cleaningLength1==0){
        cleaningLength1=parameters[2]*parameters[3];
    }
    if(cleaningLength1>thresholdForParallel2){
        parameters[14]=2;
        parameters[4]=cleaningLength1;
        cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
    }
    else{
        errNum = clEnqueueReadBuffer(*commandQueue, memObjects[5], CL_TRUE,parameters[10] * sizeof(myint),  (paramSaver1 - parameters[10]) * sizeof(myint), activeEdgesSeq+parameters[10],0, NULL, &eventE);
        
        clWaitForEvents(1, &eventE);
        
        
 
        
        cleanSequence(activeEdgesSeq, parameters+12, parameters+10, parameters+34,cleaningLength1,parameters[20] );
        
        clEnqueueWriteBuffer(*commandQueue, memObjects[5], CL_TRUE, parameters[10] * sizeof(myint), sizeof(myint) * (paramSaver1-parameters[10]), activeEdgesSeq+parameters[10],0,NULL, &eventE);
        clWaitForEvents(1, &eventE);
        
    }

    
    //*********************************************************
    //*************   Step 6: Check terminal conditions
    //*********************************************************
    //
    //  In this step we take a look whether a vertex from $B$
    //  became active or if there are no active edges. These
    //  would be the indications of the completion of the
    //  algorithm.

    keeperOfNumberOfPhantoms=parameters[29];
    parameters[29]=parameters[36];
    
    checkEnd=checkTerminalCondition(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters, sl,  pws, numMemObjects);
    
     parameters[29]=keeperOfNumberOfPhantoms;
 

    if (checkEnd==-1){
        //*********************************************************
        //*************   Step 7: Finalizing phantoms
        //*********************************************************
        //
        //  In this step we go once again over the sequence of
        //  phantoms and remove each one that has its time parameter
        //  equal to $0$.
        //  After that we clean the sequence of phantoms once again
        //  and make sure that we have the correct size of the
        //  sequence of phantoms.

        myint glSize1Ph=parameters[36];
        if((glSize1Ph%pws!=0)||(glSize1Ph==0)){
            glSize1Ph=(glSize1Ph/pws+1)*pws;
        }
        
        size_t globalWorkSizePh[1] = { glSize1Ph };
        size_t localWorkSizePh[1] = { pws };
        
        changeParameters(context,commandQueue, kernfls+6, memObjects, parameters, *numParameters);
        

        
        
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[6], 1, NULL,
                                      globalWorkSizePh, localWorkSizePh,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
        paramSaver2=(parameters[36]) * (parameters[31]);
        if(paramSaver2>parameters[30]){// The number of phantoms has exceeded the allocated memory.
            paramSaver2=parameters[30];
        }
        myint zeroSeq=0;
        if(paramSaver2==0){
            paramSaver2=  parameters[31];
            zeroSeq=1;
        }
        
        
        if(paramSaver2>thresholdForParallel2){
            parameters[4]=0;
            parameters[29]=0;
            cleanPhantomSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
            parameters[29]=parameters[36];
            changeParameters(context, commandQueue,kernfls+6, memObjects, parameters, *numParameters);
        }
        else{
            if(zeroSeq==0){
                errNum = clEnqueueReadBuffer(*commandQueue, memObjects[7], CL_TRUE,
                                             0,  (paramSaver2) * sizeof(myint), phantomSeq,
                                             0, NULL, &eventE);
                clWaitForEvents(1, &eventE);
                parameters[29]=0;
                myint bPhan=parameters[36];
                parameters[36]=0;
                cleanPhantomSequence(phantomSeq, parameters+30, parameters+29,parameters+36,bPhan,parameters[20],parameters[31]);
                parameters[29]=parameters[36];
                clEnqueueWriteBuffer(*commandQueue, memObjects[7], CL_TRUE, 0, sizeof(myint) * (paramSaver2), phantomSeq,0,NULL, &eventE);
                clWaitForEvents(1, &eventE);
                changeParameters(context, commandQueue,kernfls+6, memObjects, parameters, *numParameters);
            }
        }
        
        //*********************************************************
        //*************   Step 8: Finalizing vertices
        //*********************************************************
        //
        //  In this step of the program the sequence of active vertices
        //  is cleaned so it contains new active vertices and looses
        //  vertices that may cease to be active.
        //
        //
        //********* 8.1 Initial treatment of triggered vertices
        //  For each triggered vertex $Q$ we first check whether it was
        //  inactive before. If it was inactive then its label becomes
        //  equal to the label stored at the temporary storing location
        //  in the sequence of vertices. If it was active, its label
        //  remains unchanged. The phantoms were created and their labels
        //  are keeping track of the improved water quality that has reached
        //  the vertex $Q$.
        //  We may now clean the temporary storing location in the sequence
        //  of vertices so it now contains the symbol for emptiness (some
        //  pre-define negative number). %We also use this opportunity to
        //  erase the trigger indicator for the vertex.

        glSize1Ph=parameters[32]-parameters[9];
        if((glSize1Ph%pws!=0)||(glSize1Ph==0)){
            glSize1Ph=(glSize1Ph/pws+1)*pws;
        }
        
        
        globalWorkSizePh[0] =  glSize1Ph ;
        localWorkSizePh[0] =  pws ;
        
        changeParameters(context, commandQueue,kernfls+7, memObjects, parameters, *numParameters);
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[7], 1, NULL,
                                      globalWorkSizePh, localWorkSizePh,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        

        //********* 8.2 Merging triggered with active vertices
        //  Triggered vertices are now merged to the sequence of active vertices.
        

        
        paramSaver1=parameters[32];
        keeperOfNumberOfPhantoms=parameters[32];
        parameters[9]=0;
        parameters[32]=0;
        
        if(paramSaver1==0){
            paramSaver1=20;
        }
        
        if(paramSaver1> thresholdForParallel){
            parameters[14]=1;
            parameters[4]=paramSaver1;
            cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
            parameters[9]=parameters[32];
        }
        else{
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[4], CL_TRUE,
                                         0,  paramSaver1 * sizeof(myint), activeVerticesSeq,
                                         0, NULL, &eventE);
            
            clWaitForEvents(1, &eventE);
            cleanSequence(activeVerticesSeq, parameters+11, parameters+9, parameters+32,keeperOfNumberOfPhantoms,parameters[20]);
            parameters[9]=parameters[32];
            clEnqueueWriteBuffer(*commandQueue, memObjects[4], CL_TRUE, 0, sizeof(myint) * paramSaver1, activeVerticesSeq,0,NULL, &eventE);
            clWaitForEvents(1, &eventE);

        }
        
        
        

        //********* 8.3 Check active vertices for potential loss of activity
        //  For each active vertex $Q$ look at all edges from $Q$. If there
        //  is no active edge whose source is $Q$, then $Q$ should not be
        //  active any longer.
        
        glSize1Ph=parameters[32];
        if((glSize1Ph%pws!=0)||(glSize1Ph==0)){
            glSize1Ph=(glSize1Ph/pws+1)*pws;
        }

        globalWorkSizePh[0] =  glSize1Ph ;
        localWorkSizePh[0] =  pws ;
        
        changeParameters(context, commandQueue,kernfls+8, memObjects, parameters, *numParameters);
        
        
        
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[8], 1, NULL,
                                      globalWorkSizePh, localWorkSizePh,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        keeperOfNumberOfPhantoms=parameters[32];
        paramSaver1=parameters[32];
        if(paramSaver1==0){
            paramSaver1=20;
        }
        parameters[9]=0; parameters[32]=0;
        
        //********* 8.4 Condensing the sequence of active vertices
        //  After previous few steps some vertices may stop being active
        //  in which case they should be removed from the sequence.
        
        if(paramSaver1> thresholdForParallel){
            parameters[14]=1;
            parameters[4]=paramSaver1;

            cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
            parameters[9]=parameters[32];
        }
        else{
            
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                         0,  (*numParameters) * sizeof(myint), parameters,
                                         0, NULL, &eventE);
            
            clWaitForEvents(1, &eventE);
            paramSaver1=parameters[32];
            if(paramSaver1==0){
                paramSaver1=20;
            }
            
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[4], CL_TRUE,
                                         0,  paramSaver1 * sizeof(myint), activeVerticesSeq,
                                         0, NULL, &eventE);
            
            clWaitForEvents(1, &eventE);
            
            
            
            parameters[9]=0;
            
            condenseSequence(activeVerticesSeq,  parameters+32,parameters[20]);
            
            
            parameters[9]=parameters[32];
            changeParameters(context,commandQueue, kernfls+8, memObjects, parameters, *numParameters);
            
            clEnqueueWriteBuffer(*commandQueue, memObjects[4], CL_TRUE, 0, sizeof(myint) * paramSaver1, activeVerticesSeq,0,NULL, &eventE);
            clWaitForEvents(1, &eventE);
            
            
        }
        changeParameters(context,commandQueue, kernfls+8, memObjects, parameters, *numParameters);

        
        
        
        //*********************************************************
        //*************   Step 9: Finalizing active edges
        //*********************************************************
        //
        //  We first need to merge the triggered edges with active edges.
        
        paramSaver1=parameters[34];
        if(paramSaver1==0){
            paramSaver1=20;
        }
        
        keeperOfNumberOfPhantoms=parameters[34];
        parameters[10]=0; parameters[34]=0;
        if(keeperOfNumberOfPhantoms>thresholdForParallel2){
            parameters[14]=2;
            parameters[4]=keeperOfNumberOfPhantoms;
            cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
        }
        else{
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[5], CL_TRUE,
                                         0,  paramSaver1 * sizeof(myint), activeEdgesSeq,
                                         0, NULL, &eventE);
            clWaitForEvents(1, &eventE);
            parameters[10]=0; parameters[34]=0;
            cleanSequence(activeEdgesSeq, parameters+12, parameters+10, parameters+34,keeperOfNumberOfPhantoms,parameters[20]);
            parameters[10]=parameters[34];
            clEnqueueWriteBuffer(*commandQueue, memObjects[5], CL_TRUE, 0, sizeof(myint) * paramSaver1, activeEdgesSeq,0,NULL, &eventE);
            clWaitForEvents(1, &eventE);
        }
        

        //  Then all just used edges have to become used and their source
        //  has to be re-set so it is not equal to any of the endpoints.

        
        changeParameters(context,commandQueue,kernfls+9, memObjects, parameters, *numParameters);
        
        glSize1Ph=parameters[34];
        pws1=pws;

        if((glSize1Ph%pws1!=0)||(glSize1Ph==0)){
            glSize1Ph=(glSize1Ph/pws1+1)*pws1;
        }
        globalWorkSizePh[0] =  glSize1Ph ;
        localWorkSizePh[0] =  pws1 ;
        
        
        errNum=clEnqueueNDRangeKernel(*commandQueue, kernfls[9], 1, NULL,
                                      globalWorkSizePh, localWorkSizePh,
                                      0, NULL, &eventE);
        
        if (errNum != CL_SUCCESS){
            treatError(errNum,*context, *commandQueue, *program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
        clWaitForEvents(1, &eventE);
        
        //  Those used edges should be removed from the sequence of  active edges.
        
        errNum = clEnqueueReadBuffer(*commandQueue, memObjects[3], CL_TRUE,
                                     0,  (*numParameters) * sizeof(myint), parameters,
                                     0, NULL, &eventE);
        clWaitForEvents(1, &eventE);
        paramSaver1=parameters[34];
        if(paramSaver1==0){
            paramSaver1=20;
        }
        keeperOfNumberOfPhantoms=parameters[34];
        if(keeperOfNumberOfPhantoms>thresholdForParallel2){
            parameters[34]=0;
            parameters[14]=2;
            parameters[4]=keeperOfNumberOfPhantoms;
            cleanSequence_parallel(context, commandQueue, program, kernfls, numKernels, memObjects, parameters, numParameters,  sl,  pws, numMemObjects);
            parameters[10]=parameters[34];
        }
        else{
            errNum = clEnqueueReadBuffer(*commandQueue, memObjects[5], CL_TRUE,
                                         0,  paramSaver1 * sizeof(myint), activeEdgesSeq,
                                         0, NULL, &eventE);
            clWaitForEvents(1, &eventE);
            keeperOfNumberOfPhantoms=parameters[34];
            parameters[10]=0;
            condenseSequence(activeEdgesSeq, parameters+34,parameters[20]);
            parameters[10]=parameters[34];
            clEnqueueWriteBuffer(*commandQueue, memObjects[5], CL_TRUE, 0, sizeof(myint) * paramSaver1, activeEdgesSeq,0,NULL, &eventE);
            clWaitForEvents(1, &eventE);
        }
        changeParameters(context,commandQueue, kernfls+9, memObjects, parameters, *numParameters);
        if(parameters[10]==0){
            checkEnd=-2;
        }
    }
    return checkEnd;
}



myint shortestPathParallel(Graph* myGraph){

    myint numKernels=29;
    myint *parameters;
    myint *numberOfParameters;
    numberOfParameters=new myint;
    *numberOfParameters=61;
    myint defaultNumberForEmpty=-7;
    myint bChoice;
    parameters=new myint[*numberOfParameters];
    myint n, newDim;
    
    parameters[0]=1;parameters[1]=0;parameters[2]=5;parameters[3]=8;parameters[4]=0;parameters[5]=0;
    parameters[6]=-1;parameters[31]=5;parameters[29]=0; parameters[37]=-71; parameters[20]=defaultNumberForEmpty;
    getSetupFromFile("setup.txt", parameters+17, parameters+7, &newDim, parameters + 26, parameters+ 27, &n, parameters+39, &bChoice, parameters + 42, parameters+43);

    for(myint i=46;i<55;i++){
        parameters[i]=0;
    }
    

    myint numMemObjects=9;
    myint auxNumber=5000;

    //This auxNumber will be used to create extra memory for communication between computing cores
    //Obviously it cannont be infinity; however it must be at least as big as the number of computing cores.
    //The number appears only in percolation_parallel.cpp: createMemObjects

    string outputTemp="outputt1.txt";
    string outputAlt="example_a_out.txt";
    string inputAlt="example_a.txt";
    string informationFinish="info_fin.txt";
    ofstream fInfo;

    ostringstream st2;
    string st3;
    time_t tBeg,tEnd, pTimeKeep, npTimeKeep, pTimeKeep2;
        VtxIt begV,endV;
    cl_context context = 0;
    cl_command_queue commandQueue = 0;
    cl_program program = 0;
    cl_device_id device = 0;
    cl_kernel *kernfls;
    myint counterTime=0;
    myint indFinished=-1;

    Graph **graphP2;
    graphP2=new Graph *;
    if(parameters[45]==0){
        copyGraph(myGraph,graphP2);
    }
    kernfls=new cl_kernel[numKernels];
    myint temp;
    cl_mem *memObjects;
    memObjects= new cl_mem[numMemObjects];
    for (mysint j=0;j<numMemObjects;j++){
        memObjects[j]=0;
    }
    cl_int errNum;
    cl_event eventE;
    myint forDecrease;
    pTimeKeep=0;npTimeKeep=0;pTimeKeep2=0; tBeg=0;

    context = CreateContext(parameters[43]);
    if (context == NULL)
    {
        std::cerr << "Failed to create OpenCL context." << std::endl;
        return 1;
    }
 
    
    
    commandQueue = CreateCommandQueue(context, &device);
    if (commandQueue == NULL)
    {
        Cleanup(context, commandQueue, program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    

    program = CreateProgram(context, device, "parallel.cl");
    if (program == NULL)
    {
        Cleanup(context, commandQueue, program, kernfls, numKernels, memObjects, numMemObjects);
        return 1;
    }
    
    

    kernfls[0]=clCreateKernel(program, "terminal_condition_check", NULL);
    kernfls[1]=clCreateKernel(program, "first_step", NULL);
    kernfls[2]=clCreateKernel(program, "prep_trigger_edges_non_phantom", NULL);
    kernfls[3]=clCreateKernel(program, "prep_trigger_edges_phantom", NULL);
    kernfls[4]=clCreateKernel(program, "trigger_edges", NULL);
    kernfls[5]=clCreateKernel(program, "treatment_of_triggered_edges", NULL);
    kernfls[6]=clCreateKernel(program, "finalize_phantoms", NULL);
    kernfls[7]=clCreateKernel(program, "final_initial_treatment_of_tvs", NULL);
    kernfls[8]=clCreateKernel(program, "check_loss_of_activity_vertices", NULL);
    kernfls[9]=clCreateKernel(program, "finalize_edges", NULL);

    kernfls[10]=clCreateKernel(program, "mergeSortDecExtraStep_parallel", NULL);
    kernfls[11]=clCreateKernel(program, "mergeSortDecExtraStep_parallel_edges", NULL);
    
    kernfls[12]=clCreateKernel(program, "mergeSortDec_parallel_positionCalcs", NULL);
    kernfls[13]=clCreateKernel(program, "mergeSortDec_parallel_merging", NULL);
    kernfls[14]=clCreateKernel(program, "mergeSortDec_parallel_clean_stage1", NULL);
    kernfls[15]=clCreateKernel(program, "mergeSortDec_parallel_clean_stage2", NULL);
    kernfls[16]=clCreateKernel(program, "mergeSortDec_parallel_spot_removalSign_in_sorted", NULL);
    
    kernfls[17]=clCreateKernel(program, "mergeSortDec_parallel_positionCalcs_e", NULL);
    kernfls[18]=clCreateKernel(program, "mergeSortDec_parallel_merging_e", NULL);
    kernfls[19]=clCreateKernel(program, "mergeSortDec_parallel_clean_stage1_e", NULL);
    kernfls[20]=clCreateKernel(program, "mergeSortDec_parallel_clean_stage2_e", NULL);
    kernfls[21]=clCreateKernel(program, "mergeSortDec_parallel_spot_removalSign_in_sorted_e", NULL);
    
    
    kernfls[22]=clCreateKernel(program, "mergeSortDecBlock_parallel_positionCalcs", NULL);
    kernfls[23]=clCreateKernel(program, "mergeSortDecBlock_parallel_merging", NULL);
    kernfls[24]=clCreateKernel(program, "mergeSortDecBlockExtraStep_parallel", NULL);
    kernfls[25]=clCreateKernel(program, "mergeSortDecBlock_parallel_clean_stage1", NULL);
    kernfls[26]=clCreateKernel(program, "mergeSortDecBlock_parallel_clean_stage2", NULL);
    kernfls[27]=clCreateKernel(program, "mergeSortDecBlock_parallel_spot_removalSign_in_sorted", NULL);
    kernfls[28]=clCreateKernel(program, "mergeSort_increase_block_size_and_change_direction", NULL);

    for(myint k=0;k<numKernels;k++){
        if (kernfls[k] == NULL)
        {
            std::cerr << "Failed to create kernel " <<k<< std::endl;
            Cleanup(context, commandQueue, program, kernfls, numKernels, memObjects, numMemObjects);
            return 1;
        }
    }
    
    
    myint finished =0;
    Graph **graphP3;
    graphP3=new Graph *;
    myint currentTime, timeStep;
    myint *shortestPathSeq;
    myint shPath[1000];
    myint numTermsShP=0;
    tBeg=time(0);
    myint  counterCycl=0;
    myint algorithmChoice=parameters[42]; myint memoryCreationTime=0;
    myint * activeVerticesSeq, *activeEdgesSeq, *phantomSeq;

    while(finished==0){
        
        if(algorithmChoice==0){
            cout<<"Cycle "<<1+numTermsShP<<endl;
            copyGraph(*graphP2,graphP3);
        }
        else{
            cout<<"Creating MemObjects "<<endl;
        }
 
        
        
        tBeg=time(0);
        CreateMemObjects(context, memObjects, *graphP2, parameters,numberOfParameters, numMemObjects, auxNumber);

        for(myint i=0;i<numKernels;i++){
            setMainKernelArguments(kernfls+i, memObjects, numMemObjects);
        }

        if(algorithmChoice>0){
            tEnd=time(0);
            cout<<"Creation of memory objects is completed. It took "<<tEnd-tBeg<<" seconds. "<<endl;
            memoryCreationTime=tEnd-tBeg;
            cout<<"Writing memory objects to file."<<endl;
            cout<<"Writing completed. ";
            cout<<"Starting the timer."<<endl;
            tBeg=time(0);
        }
        
        activeVerticesSeq=new myint[parameters[11]];
        activeEdgesSeq=new myint[parameters[12]];
        phantomSeq=new myint[parameters[30]];
        
        size_t * preferred_workgroup_size, *preferred_workgroup_size02;
        
        preferred_workgroup_size=new size_t;
        preferred_workgroup_size02=new size_t;
        clGetKernelWorkGroupInfo (kernfls[0],
                                  device,
                                  CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                  sizeof(size_t),
                                  preferred_workgroup_size,
                                  NULL);
        clGetKernelWorkGroupInfo (kernfls[0], device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,  sizeof(size_t),
                                  preferred_workgroup_size02 ,NULL);
        
         size_t sl= *preferred_workgroup_size ;
        
        

        
        

        finished=0;

        indFinished=-1;
        if(algorithmChoice>0){
            counterCycl=0;
        }
        
        
        while(indFinished==-1){

            indFinished= algCycle(&context, &commandQueue, &program, kernfls, numKernels, memObjects, parameters, numberOfParameters, sl,  *preferred_workgroup_size02, numMemObjects,  activeVerticesSeq, activeEdgesSeq, phantomSeq);
            if(numTermsShP==0){counterTime++;}
 

        }

        if(indFinished!=-2){

            forDecrease=parameters[27]+parameters[8]-1;
            for(tie(begV,endV)=vertices(**graphP2);begV!=endV;++begV){
                temp=get(&VertexProperties::vCost,**graphP2,*begV);
                if(get(&VertexProperties::belongsToD,**graphP2,*begV)==2){

                    put(&VertexProperties::belongsToD,**graphP2,*begV,0);
                }
                if(get(&VertexProperties::vName,**graphP2,*begV)==indFinished){
                    if(get(&VertexProperties::belongsToD,**graphP2,*begV)==1){
                        finished=7;
                    }
                    else{
                        put(&VertexProperties::belongsToD,**graphP2,*begV,2);
                    }
                }
                if(temp>0){
                    if(temp>forDecrease){
                        temp-=forDecrease;
                    }
                    else{
                        temp=0;
                    }
                    put(&VertexProperties::vCost,**graphP2,*begV,temp);
                }
            }

            shPath[numTermsShP]=indFinished;
            numTermsShP++;
            
            parameters[0]=1;parameters[1]=0;parameters[2]=5;parameters[3]=8;parameters[4]=0;parameters[5]=0;
            parameters[6]=-1;parameters[31]=5;parameters[29]=0; parameters[37]=-71;
        }
        else{
            finished=5;
        }
        if(algorithmChoice>0){finished=7;}
        if(algorithmChoice==0){delete *graphP3;}
    
        delete[] phantomSeq;
        delete[] activeVerticesSeq;
        delete[] activeEdgesSeq;
        
    }
    
    tEnd=time(0);
    if(algorithmChoice==0){
        for(myint i=numTermsShP;i>0;i--){
            cout<<shPath[i-1]<<" ";
        }
        cout<<endl;
    }
    
    if(algorithmChoice>0){
        fInfo.open(informationFinish);
        fInfo<<"Finished. The second to last vertex on the path is "<<indFinished<<"."<<endl;
        fInfo<<"The shortest passage time: "<<counterTime<<endl;
        fInfo<<"Runtime: "<<memoryCreationTime<<" + "<< tEnd-tBeg<<endl;
        cout<<"Finished. The second to last vertex on the path is "<<indFinished<<"."<<endl;
        cout<<"The shortest passage time: "<<counterTime<<endl;
        cout<<"Runtime: "<<memoryCreationTime<<" + "<< tEnd-tBeg<<endl;
        fInfo.close();
    }

    
 

    


    Cleanup(context, commandQueue, program, kernfls, numKernels, memObjects, numMemObjects);
    delete[] memObjects;
    delete[] parameters;
    delete[] kernfls;
    delete numberOfParameters;
    delete *graphP2;
    delete graphP2;
    delete graphP3;
    return 1;
    
    
}



