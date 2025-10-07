/**********************************************************************************
* Author: Christoph Pflaum, Riccarda Scherner-Griesshammer
 *                Department Informatik Lehrstuhl 10 - Systemsimulation
 *                Friedrich-Alexander Universität Erlangen-Nürnberg
 *
*********************************************/

#include "vector.h"
#include "multilevelvector.h"
#include "../tests/testing.h"

VectorSparseG::VectorSparseG(AdaptiveSparseGrid& sg) {
    sparseGrid =  &sg;
    constructed = true;

    //number = sparseGrid->getNewDataNumber();
    length = int(sparseGrid->getMaximalOccupiedSecondTable());

    dataTableVector = new double[length];
    for(unsigned long j=0; j<length;j++)
        dataTableVector[j]=0.0;

};

VectorSparseG::VectorSparseG(AdaptiveSparseGrid* sg){
    sparseGrid = sg;
    constructed = true;

    //number = sparseGrid->getNewDataNumber();
    length = int(sparseGrid->getMaximalOccupiedSecondTable());



    dataTableVector = new double[length];
    for(unsigned long j=0; j<length;j++)
        dataTableVector[j]=0.0;


};


VectorSparseG::VectorSparseG(VectorSparseG& u) {
  constructed = true;
  sparseGrid = u.getSparseGrid();
  //number = sparseGrid->getNewDataNumber();
   length = int(sparseGrid->getMaximalOccupiedSecondTable());
   dataTableVector= new double [length];
  
  for(int i=0; i < sparseGrid->getMaximalOccupiedSecondTable(); i++){
      dataTableVector[i]=u.dataTableVector[i];
  }

};


VectorSparseG::~VectorSparseG() {
    if(constructed)
    delete[] dataTableVector;

};


void VectorSparseG::operator=(const VectorSparseG &v) {
    if (mpi_doit()) {
        AdaptiveSparseGrid_Base *sparseGrid_new = v.getSparseGrid();
        unsigned long endIndex = sparseGrid_new->maximalOccupiedSecondTable;
        //unsigned long endIndex    = sparseGrid->maximalOccupiedSecondTable;
        //double* dataTable         = sparseGrid->dataTable;
        dataInteger *secondTable = sparseGrid_new->secondTable;
        //unsigned int numberOfData = sparseGrid->numberOfData;

        //IndexDimension Idummy;
        for (unsigned long i = 0; i < endIndex; ++i) {
            if (secondTable[i] != 0) {
                IndexDimension I = sparseGrid_new->getIndexOfTable(i);
                unsigned long int k;
                if (sparseGrid->occupied(k, I)) {
                    dataTableVector[k] = v.dataTableVector[i];
                }
            }
        }
    }
}

void VectorSparseG::operator= (VectorSparseG& v) {


    if(constructed) {
        if (mpi_doit()) {


            AdaptiveSparseGrid_Base *sparseGrid_new = v.getSparseGrid();

            unsigned long endIndex = sparseGrid_new->getMaximalOccupiedSecondTable();

            for (unsigned long i = 0; i < endIndex; ++i) {
                IndexDimension I = sparseGrid_new->getIndexOfTable(i);

                if (sparseGrid_new->workonindex(i)) {
                    if(sparseGrid_new->getKey()==sparseGrid->getKey()) setValue(i,v.getValue(i));
                    else {
                        unsigned long int k;
                        if (sparseGrid->occupied(k, I)) {
                            if (workonindex(k))
                                setValue(k, v.getValue(i));
                        }
                    }

                }



            }
        }
    }
    else {

        constructed = true;

        if(sparseGrid == v.getSparseGrid()) {
            sparseGrid = v.getSparseGrid();
            //number = sparseGrid->getNewDataNumber();
            int length = int(sparseGrid->secondTableLength);
            dataTableVector = new double[length];

            for (int i = 0; i < length; i++) {
                dataTableVector[i] = v.dataTableVector[i];
            }
        }else{
            AdaptiveSparseGrid_Base* sparseGrid_new = v.getSparseGrid();
            unsigned long endIndex = sparseGrid_new-> maximalOccupiedSecondTable;
            //unsigned long endIndex    = sparseGrid->maximalOccupiedSecondTable;
            //double* dataTable         = sparseGrid->dataTable;
            dataInteger* secondTable  = sparseGrid_new->secondTable;
            //unsigned int numberOfData = sparseGrid->numberOfData;

            //IndexDimension Idummy;
            for(unsigned long i = 0;i < endIndex; ++i) {
                if (secondTable[i] != 0) {
                    if (workonindex(i)) {

                        IndexDimension I = sparseGrid_new->getIndexOfTable(i);
                        unsigned long int k;
                        if (sparseGrid->occupied(k, I)) {
                            if (workonindex(k))
                                dataTableVector[k] = v.dataTableVector[i];
                        }
                    }
                }
            }
        }


    }



}
	
bool VectorSparseG::operator==(const VectorSparseG& v)
{  
    unsigned long endIndex    = sparseGrid->maximalOccupiedSecondTable;
    dataInteger* secondTable  = sparseGrid->secondTable;
   
    for (int i= 0; i < endIndex; i++){
     
        if(secondTable[i]!=0){
        double abs = dataTableVector[i] - v.dataTableVector[i];
        
        
        
   
        if(abs<0.0) abs = -abs;
        

        if (abs < 1e-20 ){
           
            return false;          
        } 
    }
        
    }
    return true;
    
}
	

	


void VectorSparseG::Print(int level)
{
    myAssert(DimensionSparseGrid <= 2);
    int k= DimensionSparseGrid;
    if(k ==1){
        PrintIntOneD(level);
    }
    
    if (k ==2){
        PrintIntTwoD(level);
    }
}

void VectorSparseG::PrintDouble(int level) {


    if (DimensionSparseGrid == 1) {
        PrintDoubleOneD(level);
    }

    if (DimensionSparseGrid == 2) {
        PrintDoubleTwoD(level);
    }

    if (DimensionSparseGrid > 2)
        PrintDoubleTwoD(level);
}


void VectorSparseG::PrintDoubleTwoD(int level, Depth T) {
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << " integer values on grid " << endl;

    cout << "------------------------- " << endl;
    cout << endl;
    IndexDimension I;
    I = IndexOneD(startUnitInterval);
    unsigned long int indexOfData = 0;
    ////unsigned int numberOfData = sparseGrid->numberOfData;
    ////double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Depth Tlocal(I);
            if (sparseGrid->occupied(indexOfData, I) && Tlocal <= T) {
                // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";


                cout << dataTableVector[indexOfData] << "  ";
            } else {
//	      cout << " - ";
                cout << "   ";
            }
            I = I.nextRight(0, level);
        }
        cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
        I = I.nextRight(1, level);

    }
    cout << endl;
    cout << "------------------------- " << endl;
    cout << endl;

}


void VectorSparseG::PrintIntTwoD(int level) {
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << " integer values on grid " << endl;

    cout << "------------------------- " << endl;
    cout << endl;
    IndexDimension I;
    I = IndexOneD(startUnitInterval);
    unsigned long int indexOfData = 0;
    ////unsigned int numberOfData = sparseGrid->numberOfData;
    ////double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (sparseGrid->occupied(indexOfData, I)) {
                // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";	   
                cout << (int) dataTableVector[indexOfData] << "  ";
            } else {
//	      cout << " - ";
                cout << "   ";
            }
            I = I.nextRight(0, level);
        }
        cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
    I = I.nextRight(1,level);  
    
}
cout << endl;
cout << "------------------------- " << endl;
cout << endl;
    
}



	
void VectorSparseG::PrintDoubleTwoD(int level) {
    int N = 1;
    for (int i = 0; i < level; ++i) N = N * 2;
    N = N + 1;

    cout << " integer values on grid " << endl;

    cout << "------------------------- " << endl;
    cout << endl;
    IndexDimension I;
    I = IndexOneD(startUnitInterval);


    for (int d = 2; d < DimensionSparseGrid; d++)
        I.replace(d, centerUnitInterval);


    unsigned long int indexOfData = 0;
    ////unsigned int numberOfData = sparseGrid->numberOfData;
    ////double* dataTable        = sparseGrid->dataTable;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (sparseGrid->occupied(indexOfData, I)) {
                if(sparseGrid->getActiveTable()[indexOfData]){
                    double c = dataTableVector[indexOfData];
                    std::cout << c << "  ";
                }else{
                    cout << "   ";
                }
	        }
	        else {
//	      cout << " - ";
	      cout << "   ";	      
	        }
	   I = I.nextRight(0,level);
       }
        cout << endl;
        I.replace(0, IndexOneD(startUnitInterval));
    I = I.nextRight(1,level);  
    
}
cout << endl;
cout << "------------------------- " << endl;
cout << endl;
    
}


void VectorSparseG::PrintIntOneD(int level) {
   int N = 1;
   for(int i=0;i<level;++i) N = N * 2;
   N = N + 1;
   
   cout << " integer values on grid " << endl;
   cout << "------------------------- " << endl;
   IndexDimension I;
   I = IndexOneD(startUnitInterval);
   unsigned long indexOfData;
 
   //unsigned int numberOfData = sparseGrid->numberOfData;
   //double* dataTable        = sparseGrid->dataTable;
   //for(int i=0;i<N;++i) {
   int val = 0;

   for(int j=0;j<N;j++) {
       
	   if(sparseGrid->occupied(indexOfData,I)) {
	     // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";	
           
           cout << (int) dataTableVector[indexOfData] << " ";
	   }
	   else {
//	      cout << " - ";
	      cout << "   ";	      
	   }
	   I = I.nextRight(0,level);
       }
       //cout << endl; 
       //I.replace(0,IndexOneD(startUnitInterval));
       //I = I.nextRight(1,level);  
       cout << "     " << "\n";
        cout << endl;
        cout << "\n";
//}
   }
   void VectorSparseG::PrintDoubleOneD(int level) {
   int N = 1;
   for(int i=0;i<level;++i) N = N * 2;
   N = N + 1;
   
   cout << "double values on grid " << endl;
   cout << "------------------------- " << endl;
   IndexDimension I;
   I = IndexOneD(startUnitInterval);
   unsigned long indexOfData;
   //unsigned int numberOfData = sparseGrid->numberOfData;
   //double* dataTable        = sparseGrid->dataTable;
   //for(int i=0;i<N;++i) {
       for(int j=0;j<N;++j) {
	   if(sparseGrid->occupied(indexOfData,I)&& sparseGrid->workonindex(indexOfData)) {
	     // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";	      
           cout <<(double) dataTableVector[indexOfData] << " ";
	   }
	   else {
//	      cout << " - ";
	      cout << "   ";	      
	   }
	   I = I.nextRight(0,level);
       }
       //cout << endl; 
       //I.replace(0,IndexOneD(startUnitInterval));
       //I = I.nextRight(1,level);  
        cout << endl;
//}
   }
  

void VectorSparseG::PrintIndexTwoD(int d, int level) {
   int N = 1;
   for(int i=0;i<level;++i) N = N * 2;
   N = N + 1;
   
   cout << " integer values on grid " << endl;

   cout << "------------------------- " << endl;
   cout << endl;
   IndexDimension I;
   I = IndexOneD(startUnitInterval);
   unsigned long indexOfData;
   ////unsigned int numberOfData = sparseGrid->numberOfData;
   ////double* dataTable        = sparseGrid->dataTable;
  for(int i=0;i<N;++i) {
       for(int j=0;j<N;++j) {
	   if(sparseGrid->occupied(indexOfData,I)) {
	     // cout << " " << (int)dataTable[indexOfData * numberOfData + number] << " ";
//	      cout << " " << indexOfData << " ";	   
           cout << I.getIndex(d) << "  ";
	   }
	   else {
//	      cout << " - ";
	      cout << "   ";	      
	   }
	   I = I.nextRight(0,level);
       }
    cout << endl; 
    I.replace(0,IndexOneD(startUnitInterval));
    I = I.nextRight(1,level);  
    
}
cout << endl;
cout << "------------------------- " << endl;
cout << endl;
    
}


/**
 *  Prints values of vector in vtk file.
 * 
 * 
 * **/

void VectorSparseG::Print_vtk(std::ostream& Datei) {
  // Teil 0: Write information
  Datei << "# vtk DataFile Version 2.0\n"
        << "SparseGrid" << endl
        << "ASCII\n"
        << "DATASET UNSTRUCTURED_GRID\n";

  // Teil 1: Kopfzeile schreiben
  int num_total = 0;
    for(unsigned long  i=0;i<sparseGrid->maximalOccupiedSecondTable;++i) {
      if (sparseGrid->secondTable[i]!=0) num_total++;
  }
  Datei << "POINTS " << num_total << " float\n";

  
  // Teil 2: Koordinaten der Punkte ausgeben
  
  for(unsigned long  i=0;i<sparseGrid->maximalOccupiedSecondTable;++i) {
      if(sparseGrid->secondTable[i]!=0) {
	 for(int d=0;d<DimensionSparseGrid;++d) {
	     Datei <<  IndexOneD(sparseGrid->indicesSecondTable[d + i *DimensionSparseGrid]).coordinate() << "  ";
	 }
	 Datei << endl;
      }
  }
  // Teil 3: Zellen ausgeben
  Datei << "\nCELLS " << num_total << " " <<  num_total*2  << "\n";
  for(int i = 0; i < num_total; i++)
      Datei << 1 << " " << i << endl;
  Datei << "\nCELL_TYPES " << num_total << "\n";
 for(int i = 0; i < num_total; i++)
      Datei << "1" << endl;
 // Teil 4:
 Datei << "\nPOINT_DATA " << num_total << "\n";
 Datei << "SCALARS " << "SparseGrid" << " float 1\n";
 Datei << "LOOKUP_TABLE default\r\n";
 
for(unsigned long  i=0;i<sparseGrid->maximalOccupiedSecondTable;++i) {
      if (sparseGrid->secondTable[i]!=0)  Datei << dataTableVector[i] << endl;
  }
}


/**
 * Prints values of vector in gnu file. Usage for 1D or 2D grid.
 * 
 * **/


void VectorSparseG::Print_gnu(string name){
    std::ofstream Datei;

    Datei.open(name, std::ios::out);

  Datei << "# Print coordinates of adaptive sparse grid: " << endl;
  for(unsigned long  i=0;i<sparseGrid->maximalOccupiedSecondTable;++i) {
      //if(sparseGrid->getActiveTable()[i])
          if (sparseGrid->secondTable[i] != 0) {
              for (int d = 0; d < DimensionSparseGrid; ++d) {
                  Datei << IndexOneD(sparseGrid->indicesSecondTable[d + i * DimensionSparseGrid]).coordinate() << "  ";
              }
              Datei << dataTableVector[i] << endl;
              Datei << endl;
          }

  }

  Datei.close();
}






//////////////////////////////////////////////////////////
// MPI Functions
//////////////////////////////////////////////////////////







void VectorSparseG::Broadcast(int rank)
{
    
    
 
  int length = sparseGrid->secondTableLength;
  MPI_Bcast(dataTableVector,length,MPI_DOUBLE,rank,MPI_COMM_WORLD);
}




bool VectorSparseG::mpi_doit()
{
    bool doit = false;
    if(sparseGrid->mpi->all == true)doit = true;
    else {
        if((sparseGrid->mpi->myrank) == sparseGrid->mpi->rank) doit=true;
    }
    return doit;
}
    

    
void VectorSparseG::sendTo(int torank)
{
    int length = sparseGrid->secondTableLength;
    
    if (sparseGrid->mpi->rank == sparseGrid->mpi->myrank){
        
        MPI_Send(dataTableVector,length,MPI_DOUBLE,torank,0, MPI_COMM_WORLD);
      
    }else if (torank == sparseGrid->mpi->myrank){
      
        
        MPI_Status status;
        MPI_Recv(dataTableVector, length, MPI_DOUBLE,  sparseGrid->mpi->rank ,0, MPI_COMM_WORLD, &status);}

}


void VectorSparseG::ReduceSum(int rank) {
    int length = sparseGrid->secondTableLength;

    double *data = new double[length];
    MPI_Reduce(dataTableVector, data, length, MPI_DOUBLE, MPI_SUM, rank, MPI_COMM_WORLD);
    dataTableVector = data;

}

void VectorSparseG::setMultiLevelValues2(MultiLevelVector &a, Depth &T) {



    auto iter = GetNextSmallerDepthsIterator(T);
    do{
        Depth Tlocal = *iter;
        SingleDepthHashGrid& depthGrid = sparseGrid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
        const auto& mapping = depthGrid._mapPosToGridPos;
        // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
        if(depthGrid.getNumberOfEntries()>0) {
            for (size_t i = 0; i < mapping.size(); i++) {
                unsigned long k;
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                if (a.getSparseGrid()->occupied(k, I, T))
                    dataTableVector[mapping[i]] = a.getValue(k);
            }
        }
        // cout <<endl;

    }while(iter.next());

/*
 for (unsigned long i = 0; i < sparseGrid->getMaximalOccupiedSecondTable(); i++) {
        IndexDimension I = sparseGrid->getIndexOfTable(i);
        Depth Tlocal(I);
        if (Tlocal <= T) {
            unsigned long k;
            if (a.getSparseGrid()->occupied(k, I, T))
                //if(workonindex(i))
                dataTableVector[i] = a.getValue(k);
        }
    }
*/


}


void VectorSparseG::addMultiLevelValues2(MultiLevelVector &a, Depth &T) {



    auto iter = GetNextSmallerDepthsIterator(T);
    do{
        Depth Tlocal = *iter;
        SingleDepthHashGrid& depthGrid = sparseGrid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
        const auto& mapping = depthGrid._mapPosToGridPos;
        // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
        if(depthGrid.getNumberOfEntries()>0) {
            for (size_t i = 0; i < mapping.size(); i++) {
                unsigned long k;
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                if (a.getSparseGrid()->occupied(k, I, T))
                    dataTableVector[mapping[i]] += a.getValue(k);
            }
        }
        // cout <<endl;

    }while(iter.next());

/*
 for (unsigned long i = 0; i < sparseGrid->getMaximalOccupiedSecondTable(); i++) {
        IndexDimension I = sparseGrid->getIndexOfTable(i);
        Depth Tlocal(I);
        if (Tlocal <= T) {
            unsigned long k;
            if (a.getSparseGrid()->occupied(k, I, T))
                //if(workonindex(i))
                dataTableVector[i] = a.getValue(k);
        }
    }
*/


}


void VectorSparseG::setMultiLevelValues(MultiLevelVector &a, Depth &T) {

    auto iter = GetNextSmallerDepthsIterator(T);
    do{
        Depth Tlocal = *iter;
        SingleDepthHashGrid& depthGrid = sparseGrid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
        const auto& mapping = depthGrid._mapPosToGridPos;

        for (size_t i = 0; i < mapping.size(); i++)
        {
            unsigned long k;
            if(sparseGrid->workonindex(mapping[i])) {
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                if (a.getSparseGrid()->occupied(k, I, T))
                    dataTableVector[mapping[i]] = a.getValue(k);
            }
        }


    }while(iter.next());

}