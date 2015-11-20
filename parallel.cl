
//*******************************************************************************************************
//******************                                                                                *****
//******************  A parallel algorithm for constrained shortest path problem                    *****
//******************                                                                                *****
//******************  Created by Ivan Matic                                                         *****
//******************  http://www.imomath.com/maticivan                                              *****
//******************  License: GNU General Public License, version 3 (GPL-3.0)                      *****
//******************                                                                                *****
//******************  The Software is provided "as is," with all faults, defects and errors,        *****
//******************  and without warranty of any kind.                                             *****
//******************  Licensor does not warrant that the Software will be free of bugs, errors,     *****
//******************  viruses or other defects, and Licensor shall have no liability of any         *****
//******************  kind for the use of or inability to use the software, the software content    *****
//******************  or any associated service.                                                    *****
//*******************************************************************************************************


#define myint int
#define mysint int



__kernel void terminal_condition_check(__global  myint *vertices,
                             __global myint *edges,
                             __global myint *numPar,
                             __global myint *parameters,
                             __global myint *activeVertices,
                             __global myint *activeEdges,
                            __global myint *auxiliarySequence,
                            __global myint *phantoms,
                             __global myint *elementsOfB)
{
    myint gid = get_global_id(0);
    myint activeInB=0;    
    myint j=0;
    myint numElB=parameters[28];
    myint currentActivityIndex,currentPhantomIndex;
    myint sourceOfTheFinalPassage=-1;
    myint m=0;
    myint firstEdgeNumber, lastEdgeNumber;
    myint costOfTheLastEdge=-1;

    if(gid==0){
        
         while ((activeInB==0)&&(j<numElB)){
            currentActivityIndex=elementsOfB[j]+2;
            if((vertices[currentActivityIndex] % 2==1)||(vertices[currentActivityIndex]>=32)){
                activeInB=1;
                parameters[27]=vertices[elementsOfB[j]+1];
                if(vertices[elementsOfB[j]+4]>parameters[27]){
                    parameters[27]=vertices[elementsOfB[j]+4];
                }
                firstEdgeNumber=vertices[elementsOfB[j]+3];
                lastEdgeNumber=2 * parameters[1] * parameters[3];

                if(elementsOfB[j]<parameters[2]*(parameters[0]-1)){
                    lastEdgeNumber=vertices[elementsOfB[j]+3+parameters[2]];
                }
                while((sourceOfTheFinalPassage==-1)&&(firstEdgeNumber<lastEdgeNumber)){
                    if(edges[firstEdgeNumber+6]==2){
                        sourceOfTheFinalPassage =edges[firstEdgeNumber+1];
                        if(vertices[sourceOfTheFinalPassage*(parameters[2])+1]-edges[firstEdgeNumber+3]>=parameters[27]){
                            costOfTheLastEdge=edges[firstEdgeNumber+3];
                        }
                        else{
                            sourceOfTheFinalPassage=-1;
                        }
                    }
                    firstEdgeNumber+=parameters[3];
                }
                
                if(sourceOfTheFinalPassage==-1){
                    m=0;
                    while((m <parameters[29])&&(sourceOfTheFinalPassage==-1)){

                        currentPhantomIndex=m* (parameters[31]);
                        if( (phantoms[currentPhantomIndex+1]) * parameters[2]==elementsOfB[j]){

                            if(phantoms[currentPhantomIndex+3]==0){
                                sourceOfTheFinalPassage=phantoms[currentPhantomIndex];
                                costOfTheLastEdge=phantoms[currentPhantomIndex+4];
                            }
                        }
                        m++;
                    }
                }
                if(sourceOfTheFinalPassage==-1){
                    sourceOfTheFinalPassage= 0 - elementsOfB[j];
                }
                parameters[6]=sourceOfTheFinalPassage;
                parameters[8]=costOfTheLastEdge;
            }
            j++;
        }
        
    
    }
    
}


__kernel void first_step(__global  myint *vertices,
                         __global myint *edges,
                         __global myint *numPar,
                         __global myint *parameters,
                         __global myint *activeVertices,
                         __global myint *activeEdges,
                         __global myint *auxiliarySequence,
                         __global myint *phantoms,
                         __global myint *elementsOfB)
{
    myint gid = get_global_id(0);
    myint numberOfIntegersForEachVertex=parameters[2];
    myint readingLocation;
    readingLocation=gid+parameters[33];
    myint edgeLocation;
    
    myint triggeredPos=parameters[32]+gid;
    myint startPoint, sourcePoint, sourceLabel, pointLabel;
    if(readingLocation<parameters[10]){
        edgeLocation=activeEdges[readingLocation];
        edges[edgeLocation+6]=4;
        edges[edgeLocation+2]--;
        if(edges[edgeLocation+2]==0){
            edges[edgeLocation+6]=2;
            startPoint=edges[edgeLocation];
            sourcePoint=edges[edgeLocation+5];
            if(startPoint!=sourcePoint){
                //We will check whether the startPoint is a triggered vertex
                
                //We need first to find the label of the source (i.e. endPoint)
                sourceLabel=vertices[sourcePoint * numberOfIntegersForEachVertex+1];
                pointLabel=vertices[startPoint*numberOfIntegersForEachVertex +1];
                if(sourceLabel-edges[edgeLocation+3]>pointLabel){
                    activeVertices[triggeredPos]=startPoint*numberOfIntegersForEachVertex;
                }
            }
            
        }
        
    }
    
}




__kernel void prep_trigger_edges_non_phantom(__global  myint *vertices,
                                             __global myint *edges,
                                             __global myint *numPar,
                                             __global myint *parameters,
                                             __global myint *activeVertices,
                                             __global myint *activeEdges,
                                             __global myint *auxiliarySequence,
                                             __global myint *phantoms,
                                             __global myint *elementsOfB)
{
    myint gid = get_global_id(0);
    myint numberOfIntegersForEachVertex=parameters[2];
    myint readingLocation;
    readingLocation=gid+parameters[35]+parameters[9];
    myint  vertexPosition;
    myint firstEdge,  lastEdge, runningEdge, runningStatus;
    myint currentEdgeSourcePosition,currentBestLabel,relabelAttempt;
    
    if(readingLocation<parameters[32]){
        vertexPosition=activeVertices[readingLocation];
        
        firstEdge=vertices[vertexPosition+3];
        
        
        
        if(vertexPosition+3+numberOfIntegersForEachVertex<parameters[0]* numberOfIntegersForEachVertex){
            lastEdge=vertices[vertexPosition+3+numberOfIntegersForEachVertex];
        }
        else{
            lastEdge=2 * parameters[1] * parameters[3];
        }
        
        runningEdge=firstEdge;
        currentBestLabel=0;
        while(runningEdge<lastEdge){
            currentEdgeSourcePosition=(edges[runningEdge+5])*numberOfIntegersForEachVertex;
            if((edges[runningEdge+6]==2)&&(currentEdgeSourcePosition!=vertexPosition)){
                relabelAttempt=vertices[currentEdgeSourcePosition+1];
                relabelAttempt-=edges[runningEdge+3];
                if(relabelAttempt>currentBestLabel){
                    currentBestLabel=relabelAttempt;
                }
            }
            runningEdge+=parameters[3];
        }
        if(currentBestLabel<=vertices[vertexPosition+1]){
            activeVertices[readingLocation]=parameters[20];
        }
        else{
            runningStatus=vertices[vertexPosition+2];
            if(runningStatus<32){
                //The vertex was not triggered before and will be now
                runningStatus+=32;
            }
            vertices[vertexPosition+2]=runningStatus;
            // We will now look look at all edges that go through the given point again
            // Among all edges that are just used we will select one that results in the best relabeling
            // We will see which edges should be triggered and place them in the trigger list
            
            vertices[vertexPosition+4]=currentBestLabel;
            
            
        }
        
    }
    
}



__kernel void trigger_edges(__global  myint *vertices,
                                             __global myint *edges,
                                             __global myint *numPar,
                                             __global myint *parameters,
                                             __global myint *activeVertices,
                                             __global myint *activeEdges,
                                             __global myint *auxiliarySequence,
                                             __global myint *phantoms,
                                             __global myint *elementsOfB)
{
    myint gid = get_global_id(0);

    myint numberOfIntegersForEachVertex=parameters[2];
    myint numberOfIntegersForEachEdge=parameters[3];
    
    myint writingLocation;
    myint readingLocation;
    readingLocation=gid+parameters[35]+parameters[9];
    myint triggeredStart=parameters[34];
    myint endPoint, vertexPosition;
    myint firstEdge, lastEdge, runningEdge;
    myint currentLabel,relabelAttempt, alternatePosition,counter,destinationLabel;
    writingLocation=triggeredStart+2*(parameters[7]) *gid;
    if(readingLocation<parameters[32]){
        vertexPosition=activeVertices[readingLocation];
        
        firstEdge=vertices[vertexPosition+3];
        
        
        
        if(vertexPosition+numberOfIntegersForEachVertex+3<parameters[0]*numberOfIntegersForEachVertex){
            lastEdge=vertices[vertexPosition+3+numberOfIntegersForEachVertex];
        }
        else{
            lastEdge=2* parameters[1] * numberOfIntegersForEachEdge;
        }
        
        runningEdge=firstEdge;
        currentLabel=vertices[vertexPosition+4];counter=0;
        while(runningEdge<lastEdge){

            endPoint=edges[runningEdge+1];
            relabelAttempt=currentLabel-edges[runningEdge+3];
            destinationLabel=vertices[(endPoint*numberOfIntegersForEachVertex)+1];
            if(vertices[(endPoint*numberOfIntegersForEachVertex)+4]>destinationLabel){
                destinationLabel=vertices[(endPoint*numberOfIntegersForEachVertex)+4];
            }
            if(relabelAttempt>destinationLabel){
                // The edge should be triggered. We want to make sure it appears only once in the list
                // of triggered edges, and we will do that by using its first appearance in the edges sequence
                alternatePosition=runningEdge;
                if(edges[runningEdge+7]<alternatePosition){
                    alternatePosition=edges[runningEdge+7];
                }
                activeEdges[writingLocation+counter]=alternatePosition;
                counter++;
                
            }
            
            runningEdge+=numberOfIntegersForEachEdge;

        }

        
        
        
    }
    
}




__kernel void prep_trigger_edges_phantom(__global  myint *vertices,
                                  __global myint *edges,
                                  __global myint *numPar,
                                  __global myint *parameters,
                                  __global myint *activeVertices,
                                  __global myint *activeEdges,
                                  __global myint *auxiliarySequence,
                                  __global myint *phantoms,
                                  __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    myint nIntPh=parameters[31];
    myint posPh=parameters[33]+gid;
    myint timePh,readPos,destName,relAtt;
    myint triggerPos=parameters[32]+gid;
    myint vertexPos, runningStatus;
    if(posPh<parameters[29]){
        readPos=posPh* nIntPh;
        timePh=phantoms[readPos+3];
        if(timePh>0){
            timePh--;
            phantoms[readPos+3]=timePh;
        }
        if(phantoms[readPos+3]==0){
            // Edge may need to be triggered
            destName=phantoms[readPos+1];
            relAtt=phantoms[readPos+2]-phantoms[readPos+4];
            vertexPos=(parameters[2])* destName;
            if(relAtt>vertices[vertexPos+4]){
                vertices[vertexPos+4]=relAtt;
                runningStatus=vertices[vertexPos+2];
                if(runningStatus<32){
                    vertices[vertexPos+2]=runningStatus+32;
                }
                activeVertices[triggerPos]=vertexPos;
            }
            
        }
    }
    

}





__kernel void treatment_of_triggered_edges(__global  myint *vertices,
                            __global myint *edges,
                            __global myint *numPar,
                            __global myint *parameters,
                            __global myint *activeVertices,
                            __global myint *activeEdges,
                            __global myint *auxiliarySequence,
                            __global myint *phantoms,
                            __global myint *elementsOfB)
{
    myint gid = get_global_id(0);
    myint numberOfIntegersForEachVertex=parameters[2];
    myint writingLocation;
    myint readingLocation;
    readingLocation=gid+parameters[35]+parameters[10];
    myint sourceLabel;
    myint startPoint, endPoint, sourcePoint,destPoint, edgePosition, startPointLabel, endPointLabel, activityInd;
    myint edgeSource,edgeStatus,newEdgeStatus, phantomWritingPosition, alternatePosition;
    writingLocation=readingLocation+parameters[34]-parameters[10];
    
    
    if(readingLocation<parameters[34]){
        edgePosition=activeEdges[readingLocation];
        alternatePosition=edges[edgePosition+7];
        //Determining the source
        startPoint=edges[edgePosition];
        endPoint=edges[edgePosition+1];
        startPointLabel=vertices[startPoint* numberOfIntegersForEachVertex+1];
        if(startPointLabel<vertices[startPoint* numberOfIntegersForEachVertex+4]){
            startPointLabel=vertices[startPoint* numberOfIntegersForEachVertex+4];
        }
        endPointLabel=vertices[endPoint* numberOfIntegersForEachVertex+1];
        if(endPointLabel<vertices[endPoint* numberOfIntegersForEachVertex+4]){
            endPointLabel= vertices[endPoint* numberOfIntegersForEachVertex+4];
        }
        sourcePoint=startPoint;
        destPoint=endPoint;
        sourceLabel=startPointLabel;
        if(endPointLabel>startPointLabel){
            sourcePoint=endPoint;
            destPoint=startPoint;
            sourceLabel=endPointLabel;
        }
        edgeSource=edges[edgePosition+5];
        edgeStatus=edges[edgePosition+6];
        activityInd=( vertices[sourcePoint * numberOfIntegersForEachVertex+2])%2;
        newEdgeStatus=17;
        if(activityInd==0){
            newEdgeStatus=4;
            edges[edgePosition+2]=edges[edgePosition+4];
            edges[alternatePosition+2]=edges[edgePosition+4];
            if(edgeStatus==4){
                activeEdges[readingLocation]=parameters[20];
            }
            else{
                activeEdges[writingLocation]=alternatePosition;
            }
        }
        else{
            phantomWritingPosition=(parameters[36]+gid)*parameters[31];
            phantoms[phantomWritingPosition]=sourcePoint;
            phantoms[phantomWritingPosition+1]=destPoint;
            phantoms[phantomWritingPosition+2]=sourceLabel;
            phantoms[phantomWritingPosition+3]=edges[edgePosition+4];
            phantoms[phantomWritingPosition+4]=edges[edgePosition+3];

                activeEdges[readingLocation]=parameters[20];
 
        }
        

        if(newEdgeStatus!=17){
            edges[edgePosition+6]=newEdgeStatus;
            edges[alternatePosition+6]=newEdgeStatus;
            edges[edgePosition+5]=sourcePoint;
            edges[alternatePosition+5]=sourcePoint;
        }
        

    }
    
}



__kernel void finalize_phantoms(__global  myint *vertices,
                                __global myint *edges,
                                __global myint *numPar,
                                __global myint *parameters,
                                __global myint *activeVertices,
                                __global myint *activeEdges,
                                __global myint *auxiliarySequence,
                                __global myint *phantoms,
                                __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    myint nIntPh=parameters[31];
    myint posPh=gid * nIntPh;
    myint remSi=parameters[20];
    myint instIf,negInstIf;
    myint firstIf;
    firstIf=(gid<parameters[36]);
    myint farAway=parameters[30]-2*nIntPh;
    posPh=firstIf * posPh + (1-firstIf)*farAway;
    

    instIf=(phantoms[posPh+3]==0);
    negInstIf=1-instIf;
    phantoms[posPh]=instIf*remSi + negInstIf*phantoms[posPh];
    phantoms[posPh+1]=instIf*remSi + negInstIf*phantoms[posPh+1];
    phantoms[posPh+2]=instIf*remSi + negInstIf*phantoms[posPh+2];
    phantoms[posPh+3]=instIf*remSi + negInstIf*phantoms[posPh+3];
    phantoms[posPh+4]=instIf*remSi + negInstIf*phantoms[posPh+4];

    
}






__kernel void final_initial_treatment_of_tvs(__global  myint *vertices,
                                             __global myint *edges,
                                             __global myint *numPar,
                                             __global myint *parameters,
                                             __global myint *activeVertices,
                                             __global myint *activeEdges,
                                             __global myint *auxiliarySequence,
                                             __global myint *phantoms,
                                             __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    myint lastV=parameters[32];

    myint activityInd, readPos;
    myint vReadPos;
    readPos=parameters[9]+gid;
    if(readPos<lastV){
        vReadPos=activeVertices[readPos];
        activityInd=vertices[vReadPos+2];
        if(activityInd % 2==0){
            if(vertices[vReadPos+1]<vertices[vReadPos+4]){
                vertices[vReadPos+1]=vertices[vReadPos+4];
            }
            activityInd+=1;
        }
        if(activityInd>=32){
            vertices[vReadPos+2]=activityInd-32;
        }
    }
    
    
}

__kernel void check_loss_of_activity_vertices(__global  myint *vertices,
                                              __global myint *edges,
                                              __global myint *numPar,
                                              __global myint *parameters,
                                              __global myint *activeVertices,
                                              __global myint *activeEdges,
                                              __global myint *auxiliarySequence,
                                              __global myint *phantoms,
                                              __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    myint nIntV=parameters[2];
    myint lastV=parameters[32];


    myint readPos, vActive;
    myint vertexPosition, firstEdge, lastEdge, runningEdge,  edgeActivity, edgeSource,activeOuts;
    readPos=gid;
    if(readPos<lastV){

        vertexPosition=activeVertices[readPos];
        
        firstEdge=vertices[vertexPosition+3];
        

        
        if(vertexPosition+nIntV+3<parameters[0]*nIntV){
            lastEdge=vertices[vertexPosition+3+nIntV];
        }
        else{
            lastEdge=2*parameters[1] * parameters[3];
        }
        
        runningEdge=firstEdge;
        activeOuts=0;
 
        while(runningEdge<lastEdge){
            edgeActivity=edges[runningEdge+6];
            edgeSource=edges[runningEdge+5];
            if(((edgeActivity==3)||(edgeActivity==4))&&(edgeSource* parameters[2]==vertexPosition)){
                activeOuts++;
            }

            runningEdge+=parameters[3];
        }
        if(activeOuts==0){
            vActive=vertices[vertexPosition+2];
            vertices[vertexPosition+2]=0;
            activeVertices[readPos]=parameters[20];
            
        }
        
    }
    
    
}


__kernel void finalize_edges(__global  myint *vertices,
                                              __global myint *edges,
                                              __global myint *numPar,
                                              __global myint *parameters,
                                              __global myint *activeVertices,
                                              __global myint *activeEdges,
                                              __global myint *auxiliarySequence,
                                              __global myint *phantoms,
                                              __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    myint lastE =parameters[34];
    myint readPos;
    myint edgePosition,   edgeActivity;
    readPos=gid;
    if(readPos<lastE){
        edgePosition=activeEdges[readPos];
        edgeActivity=edges[edgePosition+6];
        if((edgeActivity==1)||(edgeActivity==2)){
            edges[edgePosition+6]=1;
            edges[edgePosition+5]=-1;
            edges[edgePosition+2]=edges[edgePosition+4];
            activeEdges[readPos]=parameters[20];
        }
        if(edgeActivity==4){
            edges[edgePosition+6]=3;
        }
        
        
    }
    
    
}

__kernel void mergeSortDec_parallel_positionCalcs(__global  myint *vertices,
                                                  __global myint *edges,
                                                  __global myint *numPar,
                                                  __global myint *parameters,
                                                  __global myint *activeVertices,
                                                  __global myint *activeEdges,
                                                  __global myint *auxiliarySequence,
                                                  __global myint *phantoms,
                                                  __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;
    myint cleaningLength=parameters[4];
    myint blockOld=parameters[38];
    myint direction=parameters[39];
    myint blockNew= 2 * blockOld;
    cleanFrom=parameters+32;
    
    myint numberOfElementsInAux=parameters[13];
    if(direction==0){
        sequence=activeVertices + (*cleanFrom);
    }
    else{
        sequence=auxiliarySequence + (*cleanFrom);
    }
    

    myint secretStorageLocation=numberOfElementsInAux-1-gid;
    
    myint myPositionInNewBlock=gid %blockNew;
    myint firstElementOfNewBlock=(gid/blockNew)*blockNew;
    myint myPositionInMyBlock;
    myint leftEnd, rightEnd;
    myint lOR=1;
    myPositionInMyBlock=myPositionInNewBlock-blockOld;
    leftEnd=firstElementOfNewBlock;
    if(myPositionInNewBlock<blockOld){
        myPositionInMyBlock+=blockOld;
        leftEnd+=blockOld;
        lOR=0;
    }
    
    rightEnd=leftEnd+blockOld;
    if(rightEnd>=cleaningLength){
        rightEnd=cleaningLength;
    }
    if(leftEnd>=cleaningLength){
        leftEnd=cleaningLength;
    }
    myint tempRightEnd=rightEnd;
    myint absLeftEnd=leftEnd;
    myint middleTerm;
    if(gid<cleaningLength){
        myint currentElement=sequence[gid];
        if(lOR==1){
            currentElement--;
        }
        
        while(leftEnd+1<tempRightEnd){
            middleTerm=(leftEnd+tempRightEnd)/2;
            if(sequence[middleTerm]>currentElement){
                leftEnd=middleTerm;
            }
            else{
                tempRightEnd=middleTerm;
            }
        }
        
        if((leftEnd<cleaningLength)&&(sequence[leftEnd]>currentElement)){
            auxiliarySequence[secretStorageLocation]=myPositionInMyBlock+ leftEnd-absLeftEnd+1 ;

        }
        else{
            auxiliarySequence[secretStorageLocation]=myPositionInMyBlock;

        }
    }
    
    
    
    
}

__kernel void mergeSortDec_parallel_merging(__global  myint *vertices,
                                            __global myint *edges,
                                            __global myint *numPar,
                                            __global myint *parameters,
                                            __global myint *activeVertices,
                                            __global myint *activeEdges,
                                            __global myint *auxiliarySequence,
                                            __global myint *phantoms,
                                            __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;
    __global myint *auxSeq;
    myint cleaningLength=parameters[4];
    myint blockOld=parameters[38];
    myint direction=parameters[39];
    myint blockNew= 2 * blockOld;
    cleanFrom=parameters+32;
    
    myint numberOfElementsInAux=parameters[13];
    
    if(direction==0){
        sequence=activeVertices + (*cleanFrom);
        auxSeq=auxiliarySequence + (*cleanFrom);
    }
    else{
        auxSeq=activeVertices + (*cleanFrom);
        sequence=auxiliarySequence + (*cleanFrom);
    }
    
    myint secretStorageLocation=numberOfElementsInAux-1-gid;
    
    
    myint firstElementOfNewBlock=(gid/blockNew)*blockNew;
    
    if(gid<cleaningLength){
        myint currentElement=sequence[gid];
        myint shift=auxiliarySequence[secretStorageLocation];
        auxSeq[firstElementOfNewBlock+ shift]=currentElement;
        
    }
    
    
    
    
}


__kernel void mergeSortDec_parallel_clean_stage1(__global  myint *vertices,
                                                      __global myint *edges,
                                                      __global myint *numPar,
                                                      __global myint *parameters,
                                                      __global myint *activeVertices,
                                                      __global myint *activeEdges,
                                                      __global myint *auxiliarySequence,
                                                      __global myint *phantoms,
                                                      __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;
    myint cleaningLength=parameters[4];

    cleanFrom=parameters+32;
    
    sequence=activeVertices + (*cleanFrom);
    myint numberOfElementsInAux=parameters[13];
    
    
    if((gid>0)&&(gid<cleaningLength)){
        myint secretStorageLocation=numberOfElementsInAux-1-gid;
        auxiliarySequence[secretStorageLocation]=0;
        if(sequence[gid]==sequence[gid-1]){
            auxiliarySequence[secretStorageLocation]=1;
        }
        
    }
    
    
    
    
}
__kernel void mergeSortDec_parallel_clean_stage2(__global  myint *vertices,
                                                      __global myint *edges,
                                                      __global myint *numPar,
                                                      __global myint *parameters,
                                                      __global myint *activeVertices,
                                                      __global myint *activeEdges,
                                                      __global myint *auxiliarySequence,
                                                      __global myint *phantoms,
                                                      __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;

    myint cleaningLength=parameters[4];
    myint removalSign=parameters[20];

    cleanFrom=parameters+32;
    
    sequence=activeVertices + (*cleanFrom);
    myint numberOfElementsInAux=parameters[13];
    
    
    if((gid>0)&&(gid<cleaningLength)){
        myint secretStorageLocation=numberOfElementsInAux-1-gid;
        if(auxiliarySequence[secretStorageLocation]==1){
            sequence[gid]=removalSign;
        }
        
    }
    
    
    
    
}

__kernel void mergeSortDec_parallel_spot_removalSign_in_sorted(__global  myint *vertices,
                                                      __global myint *edges,
                                                      __global myint *numPar,
                                                      __global myint *parameters,
                                                      __global myint *activeVertices,
                                                      __global myint *activeEdges,
                                                      __global myint *auxiliarySequence,
                                                      __global myint *phantoms,
                                                      __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;

    __global myint *cleanFrom;
    myint cleaningLength=parameters[4];
    myint removalSign=parameters[20];

    cleanFrom=parameters+32;
    
    sequence=activeVertices + (*cleanFrom);


    myint tempBeg=0;
    myint tempEnd= cleaningLength;
    myint tempMid;
    if(gid==0){
        while(tempBeg+1<tempEnd){
            tempMid=(tempBeg+tempEnd)/2;
            if(sequence[tempMid]==removalSign){
                tempEnd=tempMid;
            }
            else{
                tempBeg=tempMid;
            }
        }
        if(sequence[0]!=removalSign){
            parameters[39]=*cleanFrom+ tempEnd;
        }
        else{
            parameters[39]=*cleanFrom;
        }
        
    }

    
    
}



__kernel void mergeSortDec_parallel_positionCalcs_e(__global  myint *vertices,
                                                    __global myint *edges,
                                                    __global myint *numPar,
                                                    __global myint *parameters,
                                                    __global myint *activeVertices,
                                                    __global myint *activeEdges,
                                                    __global myint *auxiliarySequence,
                                                    __global myint *phantoms,
                                                    __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;
    myint cleaningLength=parameters[4];
    myint blockOld=parameters[38];
    myint direction=parameters[39];
    myint blockNew= 2 * blockOld;
    
    
    cleanFrom=parameters+34;
    
    myint numberOfElementsInAux=parameters[13];
    if(direction==0){
        sequence=activeEdges + (*cleanFrom);
    }
    else{
        sequence=auxiliarySequence + (*cleanFrom);
    }
    
    myint secretStorageLocation=numberOfElementsInAux-1-gid;
    
    myint myPositionInNewBlock=gid %blockNew;
    myint firstElementOfNewBlock=(gid/blockNew)*blockNew;
    myint myPositionInMyBlock;
    myint leftEnd, rightEnd;
    myint lOR=1;
    myPositionInMyBlock=myPositionInNewBlock-blockOld;
    leftEnd=firstElementOfNewBlock;
    if(myPositionInNewBlock<blockOld){
        myPositionInMyBlock+=blockOld;
        leftEnd+=blockOld;
        lOR=0;
    }
    
    rightEnd=leftEnd+blockOld;
    if(rightEnd>=cleaningLength){
        rightEnd=cleaningLength;
    }
    if(leftEnd>=cleaningLength){
        leftEnd=cleaningLength;
    }
    myint tempRightEnd=rightEnd;
    myint absLeftEnd=leftEnd;
    myint middleTerm;
    if(gid<cleaningLength){
        myint currentElement=sequence[gid];
        if(lOR==1){
            currentElement--;
        }
        
        while(leftEnd+1<tempRightEnd){
            middleTerm=(leftEnd+tempRightEnd)/2;
            if(sequence[middleTerm]>currentElement){
                leftEnd=middleTerm;
            }
            else{
                tempRightEnd=middleTerm;
            }
        }
        
        if((leftEnd<cleaningLength)&&(sequence[leftEnd]>currentElement)){
            auxiliarySequence[secretStorageLocation]=myPositionInMyBlock+ leftEnd-absLeftEnd+1 ;
            
        }
        else{
            auxiliarySequence[secretStorageLocation]=myPositionInMyBlock;
            
        }
    }
    
    
    
    
}

__kernel void mergeSortDec_parallel_merging_e(__global  myint *vertices,
                                              __global myint *edges,
                                              __global myint *numPar,
                                              __global myint *parameters,
                                              __global myint *activeVertices,
                                              __global myint *activeEdges,
                                              __global myint *auxiliarySequence,
                                              __global myint *phantoms,
                                              __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;

    __global myint *cleanFrom;
    __global myint *auxSeq;
    myint cleaningLength=parameters[4];

    myint blockOld=parameters[38];
    myint direction=parameters[39];
    myint blockNew= 2 * blockOld;

    cleanFrom=parameters+34;
    
    myint numberOfElementsInAux=parameters[13];
    
    if(direction==0){
        sequence=activeEdges + (*cleanFrom);
        auxSeq=auxiliarySequence + (*cleanFrom);
    }
    else{
        auxSeq=activeEdges + (*cleanFrom);
        sequence=auxiliarySequence + (*cleanFrom);
    }
    
    
    myint secretStorageLocation=numberOfElementsInAux-1-gid;
    
    
    myint firstElementOfNewBlock=(gid/blockNew)*blockNew;
    
    if(gid<cleaningLength){
        myint currentElement=sequence[gid];
        myint shift=auxiliarySequence[secretStorageLocation];
        auxSeq[firstElementOfNewBlock+ shift]=currentElement;
        
    }
    
    
    
    
}


__kernel void mergeSortDec_parallel_clean_stage1_e(__global  myint *vertices,
                                                   __global myint *edges,
                                                   __global myint *numPar,
                                                   __global myint *parameters,
                                                   __global myint *activeVertices,
                                                   __global myint *activeEdges,
                                                   __global myint *auxiliarySequence,
                                                   __global myint *phantoms,
                                                   __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;

    __global myint *cleanFrom;

    myint cleaningLength=parameters[4];

    cleanFrom=parameters+34;
    
    sequence=activeEdges + (*cleanFrom);
    myint numberOfElementsInAux=parameters[13];
    
    
    if((gid>0)&&(gid<cleaningLength)){
        myint secretStorageLocation=numberOfElementsInAux-1-gid;
        auxiliarySequence[secretStorageLocation]=0;
        if(sequence[gid]==sequence[gid-1]){
            auxiliarySequence[secretStorageLocation]=1;
        }
        
    }
    
    
    
    
}
__kernel void mergeSortDec_parallel_clean_stage2_e(__global  myint *vertices,
                                                   __global myint *edges,
                                                   __global myint *numPar,
                                                   __global myint *parameters,
                                                   __global myint *activeVertices,
                                                   __global myint *activeEdges,
                                                   __global myint *auxiliarySequence,
                                                   __global myint *phantoms,
                                                   __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;

    __global myint *cleanFrom;

    myint cleaningLength=parameters[4];
    myint removalSign=parameters[20];

    cleanFrom=parameters+34;
    
    sequence=activeEdges + (*cleanFrom);
    myint numberOfElementsInAux=parameters[13];
    
    
    if((gid>0)&&(gid<cleaningLength)){
        myint secretStorageLocation=numberOfElementsInAux-1-gid;
        if(auxiliarySequence[secretStorageLocation]==1){
            sequence[gid]=removalSign;
        }
        
    }
    
    
    
    
}

__kernel void mergeSortDec_parallel_spot_removalSign_in_sorted_e(__global  myint *vertices,
                                                                 __global myint *edges,
                                                                 __global myint *numPar,
                                                                 __global myint *parameters,
                                                                 __global myint *activeVertices,
                                                                 __global myint *activeEdges,
                                                                 __global myint *auxiliarySequence,
                                                                 __global myint *phantoms,
                                                                 __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;

    __global myint *cleanFrom;

    myint cleaningLength=parameters[4];
    myint removalSign=parameters[20];

    cleanFrom=parameters+34;
    
    sequence=activeEdges + (*cleanFrom);
    
    myint tempBeg=0;
    myint tempEnd= cleaningLength;
    myint tempMid;
    if(gid==0){
        while(tempBeg+1<tempEnd){
            tempMid=(tempBeg+tempEnd)/2;
            if(sequence[tempMid]==removalSign){
                tempEnd=tempMid;
            }
            else{
                tempBeg=tempMid;
            }
        }
        if(sequence[0]!=removalSign){
            parameters[39]=*cleanFrom+ tempEnd;
        }
        else{
            parameters[39]=*cleanFrom;
        }
        
    }
    
    
    
}

















__kernel void mergeSortDecExtraStep_parallel(__global  myint *vertices,
                                             __global myint *edges,
                                             __global myint *numPar,
                                             __global myint *parameters,
                                             __global myint *activeVertices,
                                             __global myint *activeEdges,
                                             __global myint *auxiliarySequence,
                                             __global myint *phantoms,
                                             __global myint *elementsOfB)
{
    myint gid=get_global_id(0);

    if(gid<parameters[4] ){
        myint position= parameters[32]+ gid;
        activeVertices[position]=auxiliarySequence[position];
    }
    
    
    
}



__kernel void mergeSortDecExtraStep_parallel_edges(__global  myint *vertices,
                                                   __global myint *edges,
                                                   __global myint *numPar,
                                                   __global myint *parameters,
                                                   __global myint *activeVertices,
                                                   __global myint *activeEdges,
                                                   __global myint *auxiliarySequence,
                                                   __global myint *phantoms,
                                                   __global myint *elementsOfB)
{
    myint gid=get_global_id(0);

    
    if(gid< parameters[4] ){
        myint position= parameters[34] + gid;
        activeEdges[position]=auxiliarySequence[position];
    }
    
    
    
}






__kernel void mergeSortDecBlock_parallel_positionCalcs(__global  myint *vertices,
                                                  __global myint *edges,
                                                  __global myint *numPar,
                                                  __global myint *parameters,
                                                  __global myint *activeVertices,
                                                  __global myint *activeEdges,
                                                  __global myint *auxiliarySequence,
                                                  __global myint *phantoms,
                                                  __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;
    __global myint *auxSeq;
    myint cleaningLength=parameters[4];

    myint blockOld=parameters[38];
    myint direction=parameters[39];
    myint blockNew= 2 * blockOld;
    myint nThPh=parameters[31];
    

    
    

    cleanFrom=parameters+36;
    
    myint numberOfElementsInAux=parameters[13];
    if(direction==0){
        sequence=phantoms + (*cleanFrom) * nThPh;
        auxSeq=auxiliarySequence + (*cleanFrom) * nThPh;
    }
    else{
        auxSeq=phantoms + (*cleanFrom) * nThPh;
        sequence=auxiliarySequence + (*cleanFrom) * nThPh;
    }
    
    myint secretStorageLocation=numberOfElementsInAux-1-gid;
    
    myint myPositionInNewBlock=gid %blockNew;
    myint firstElementOfNewBlock=(gid/blockNew)*blockNew;
    myint myPositionInMyBlock;
    myint leftEnd, rightEnd, leftEndMult;
    myint lOR=1;
    myPositionInMyBlock=myPositionInNewBlock-blockOld;
    leftEnd=firstElementOfNewBlock;
    if(myPositionInNewBlock<blockOld){
        myPositionInMyBlock+=blockOld;
        leftEnd+=blockOld;
        lOR=0;
    }
    
    rightEnd=leftEnd+blockOld;
    if(rightEnd>=cleaningLength){
        rightEnd=cleaningLength;
    }
    if(leftEnd>=cleaningLength){
        leftEnd=cleaningLength;
    }
    myint tempRightEnd=rightEnd;
    myint absLeftEnd=leftEnd;
    myint middleTermReal, middleTerm;
    if(gid<cleaningLength){
        myint currentElement1=sequence[gid * nThPh+1];
        myint currentElement2=sequence[gid * nThPh+2];
        myint currentElement3=sequence[gid * nThPh+3];
        myint currentElement4=sequence[gid * nThPh+4];
        myint compRes;
        // compRes will be 1 if the middle Term is bigger than the current element in the phantom binary comparison relation
        if(lOR==1){
            currentElement2--;
        }
        
        while(leftEnd+1<tempRightEnd){
            middleTermReal=(leftEnd+tempRightEnd)/2;
            middleTerm=middleTermReal* nThPh;
            compRes=0;
            if(sequence[middleTerm+1]>currentElement1){
                compRes=1;
            }
            else{
                if(sequence[middleTerm+1]==currentElement1){
                    if(sequence[middleTerm+3]>currentElement3){
                        compRes=1;
                    }
                    else{
                        if(sequence[middleTerm+3]==currentElement3){
                            if(sequence[middleTerm+2]-sequence[middleTerm+4]>currentElement2-currentElement4){
                                compRes=1;
                            }
                        }
                    }
                }
            }
            if(compRes==1){
                leftEnd=middleTermReal;
            }
            else{
                tempRightEnd=middleTermReal;
            }
        }
        compRes=0;
        if(leftEnd<cleaningLength){
            leftEndMult=leftEnd*nThPh;
            if(sequence[leftEndMult+1]>currentElement1){
                compRes=1;
            }
            else{
                if(sequence[leftEndMult+1]==currentElement1){
                    if(sequence[leftEndMult+3]>currentElement3){
                        compRes=1;
                    }
                    else{
                        if(sequence[leftEndMult+3]==currentElement3){
                            if(sequence[leftEndMult+2]-sequence[leftEndMult+4]>currentElement2-currentElement4){
                                compRes=1;
                            }
                        }
                    }
                }

            }
        }
        if((leftEnd<cleaningLength)&&(compRes==1)){
            auxiliarySequence[secretStorageLocation]=myPositionInMyBlock+ leftEnd-absLeftEnd+1 ;
            
        }
        else{
            auxiliarySequence[secretStorageLocation]=myPositionInMyBlock;
            
        }
    }
    
    
    
    
}

__kernel void mergeSortDecBlock_parallel_merging(__global  myint *vertices,
                                            __global myint *edges,
                                            __global myint *numPar,
                                            __global myint *parameters,
                                            __global myint *activeVertices,
                                            __global myint *activeEdges,
                                            __global myint *auxiliarySequence,
                                            __global myint *phantoms,
                                            __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;
    __global myint *auxSeq;
    myint cleaningLength=parameters[4];

    myint blockOld=parameters[38];
    myint direction=parameters[39];
    myint blockNew= 2 * blockOld;
    myint nThPh=parameters[31];

    
    
    cleanFrom=parameters+36;
    
    myint numberOfElementsInAux=parameters[13];
    
    if(direction==0){
        sequence=phantoms + (*cleanFrom) * nThPh;
        auxSeq=auxiliarySequence + (*cleanFrom) * nThPh;
    }
    else{
        auxSeq=phantoms + (*cleanFrom) * nThPh;
        sequence=auxiliarySequence + (*cleanFrom) * nThPh;
    }
    
    myint secretStorageLocation=numberOfElementsInAux-1-gid;
    
    
    myint firstElementOfNewBlock=(gid/blockNew)*blockNew;
    
    if(gid<cleaningLength){
        myint currentElement0=sequence[gid*nThPh];
        myint currentElement1=sequence[gid*nThPh+1];
        myint currentElement2=sequence[gid*nThPh+2];
        myint currentElement3=sequence[gid*nThPh+3];
        myint currentElement4=sequence[gid*nThPh+4];
        myint shift=auxiliarySequence[secretStorageLocation];
        auxSeq[(firstElementOfNewBlock+ shift)*nThPh]=currentElement0;
        auxSeq[(firstElementOfNewBlock+ shift)*nThPh+1]=currentElement1;
        auxSeq[(firstElementOfNewBlock+ shift)*nThPh+2]=currentElement2;
        auxSeq[(firstElementOfNewBlock+ shift)*nThPh+3]=currentElement3;
        auxSeq[(firstElementOfNewBlock+ shift)*nThPh+4]=currentElement4;
        
    }
    
    
    
    
}




__kernel void mergeSortDecBlockExtraStep_parallel(__global  myint *vertices,
                                             __global myint *edges,
                                             __global myint *numPar,
                                             __global myint *parameters,
                                             __global myint *activeVertices,
                                             __global myint *activeEdges,
                                             __global myint *auxiliarySequence,
                                             __global myint *phantoms,
                                             __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;

    myint cleaningLength=parameters[4];

    myint nThPh=parameters[31];
    
    
    sequence=phantoms;
    cleanFrom=parameters+36;

    
    myint position= *cleanFrom + gid;
    if(position<*cleanFrom + cleaningLength ){
        for(myint i=0;i<nThPh;i++){
            sequence[position* nThPh+i ]=auxiliarySequence[position*nThPh+i ];
        }
    }
    
    
    
}




__kernel void mergeSortDecBlock_parallel_clean_stage1(__global  myint *vertices,
                                                 __global myint *edges,
                                                 __global myint *numPar,
                                                 __global myint *parameters,
                                                 __global myint *activeVertices,
                                                 __global myint *activeEdges,
                                                 __global myint *auxiliarySequence,
                                                 __global myint *phantoms,
                                                 __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
 
    __global myint *cleanFrom;

    myint cleaningLength=parameters[4];

    myint nThPh=parameters[31];

 
    cleanFrom=parameters+36;
    
    sequence=phantoms + (*cleanFrom)*nThPh;
    myint numberOfElementsInAux=parameters[13];
    

    if((gid>0)&&(gid<cleaningLength)){
        myint secretStorageLocation=numberOfElementsInAux-1-gid;
        auxiliarySequence[secretStorageLocation]=0;
        myint equalityIndicator=0;
        if((sequence[gid*nThPh+1]==sequence[(gid-1)*nThPh+1]) && (sequence[gid*nThPh+3]==sequence[(gid-1)*nThPh+3])){
            equalityIndicator=1;
        }
 
        auxiliarySequence[secretStorageLocation]=equalityIndicator;
        
        
    }
    
    
    
    
}
__kernel void mergeSortDecBlock_parallel_clean_stage2(__global  myint *vertices,
                                                 __global myint *edges,
                                                 __global myint *numPar,
                                                 __global myint *parameters,
                                                 __global myint *activeVertices,
                                                 __global myint *activeEdges,
                                                 __global myint *auxiliarySequence,
                                                 __global myint *phantoms,
                                                 __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;
    __global myint *cleanFrom;
    myint cleaningLength=parameters[4];
    myint removalSign=parameters[20];

    myint nThPh=parameters[31];


    cleanFrom=parameters+36;
    
    sequence=phantoms + (*cleanFrom)*nThPh;
    myint numberOfElementsInAux=parameters[13];
    
    
    if((gid>0)&&(gid<cleaningLength)){
        myint secretStorageLocation=numberOfElementsInAux-1-gid;
        if(auxiliarySequence[secretStorageLocation]==1){
            for(myint i=0;i<nThPh;i++){
                sequence[gid*nThPh+i]=removalSign;
            }
        }
        
    }
    
    
    
    
}

__kernel void mergeSortDecBlock_parallel_spot_removalSign_in_sorted(__global  myint *vertices,
                                                               __global myint *edges,
                                                               __global myint *numPar,
                                                               __global myint *parameters,
                                                               __global myint *activeVertices,
                                                               __global myint *activeEdges,
                                                               __global myint *auxiliarySequence,
                                                               __global myint *phantoms,
                                                               __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    
    __global myint *sequence;

    __global myint *cleanFrom;

    myint cleaningLength=parameters[4];
    myint removalSign=parameters[20];

    myint nThPh=parameters[31];

 
    
    

    cleanFrom=parameters+36;
    
    sequence=phantoms + (*cleanFrom)*nThPh;
    myint tempBeg=0;
    myint tempEnd= cleaningLength;
    myint tempMid;
    if(gid==0){
        while(tempBeg+1<tempEnd){
            tempMid=(tempBeg+tempEnd)/2;
            if(sequence[tempMid*nThPh+1]==removalSign){
                tempEnd=tempMid;
            }
            else{
                tempBeg=tempMid;
            }
        }
        if(sequence[1]!=removalSign){
            parameters[39]=*cleanFrom+ tempEnd;
        }
        else{
            parameters[39]=*cleanFrom;
        }
        
    }
    
    
    
}

__kernel void mergeSort_increase_block_size_and_change_direction(__global  myint *vertices,
                                                                    __global myint *edges,
                                                                    __global myint *numPar,
                                                                    __global myint *parameters,
                                                                    __global myint *activeVertices,
                                                                    __global myint *activeEdges,
                                                                    __global myint *auxiliarySequence,
                                                                    __global myint *phantoms,
                                                                    __global myint *elementsOfB)
{
    myint gid=get_global_id(0);
    if(gid==0){
        parameters[39]=parameters[39]+1;
        if(parameters[39]>1){
            parameters[39]=0;
        }
        parameters[38]= parameters[38]*2;
    }
    
    
    
}
