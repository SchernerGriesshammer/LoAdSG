//
// Created by to35jepo on 6/6/24.
//

#ifndef RUN_LOCALSTIFFNESSMATRICES_H
#define RUN_LOCALSTIFFNESSMATRICES_H


#include <queue>
#include <unordered_set>
#include "../mympi.h"
#include "../stencils/Stencil.h"
#include "../cells/CellStructure.h"


double bytesToGigabytes(unsigned long long bytes);


struct Container {
    int id;
    int totalSize;

    // This operator will allow the priority_queue to compare containers based on their totalSize.
    bool operator>(const Container &other) const {
        return totalSize > other.totalSize;
    }
    std::vector<std::pair<Depth,int>> boxes;

    bool full=false;
};



// Custom hash function for Depth
struct DepthHash {
    std::size_t operator()(const Depth& D) const {
        unsigned long value = D.at(0);
        for (int d = 1; d < DimensionSparseGrid; ++d) {
            value = value + D.at(d) * PrimeNumbers::getPrimeForHash(d);
        }
        return value;
    }
};

// Custom equality comparison function for Depth
struct DepthEqual {
    bool operator()(const Depth& lhs, const Depth& rhs) const {
        // Implement an equality comparison function for Depth here
        for(int d=0; d<DimensionSparseGrid; d++){
            if(lhs.at(d)!=rhs.at(d)) return false;
        }
        return true;
    }
};


class DistributedDepthsHashtable{
    std::unordered_map<Depth,int, DepthHash, DepthEqual> _map;

public:
    DistributedDepthsHashtable()=default;
    DistributedDepthsHashtable(Container* containers_sorted, int n);
    int getNodeForDepth(const Depth &D);
    void setNodeDepth(int i, Depth& D){_map.insert({D,i});};

    std::unordered_map<Depth,int, DepthHash, DepthEqual>* getMap(){
        return &_map;
    }

};



void distribute_LocalStiffnessMatrices(int numberofprocesses, AdaptiveSparseGrid &sg, Container* containers_sorted);






class LocalStiffnessMatrices {

protected:
    class LocalStiffnessMatrixFixedDepthSymmetric;

public:
    //ghost functions
    virtual inline void initialize(Depth &T_){};
    virtual inline void applyStencilOnCell(CellDimension& cell,VectorSparseG& input, VectorSparseG& output){};
    virtual inline double returnValue(const IndexDimension &Index, const MultiDimCompass &mc){return 0.0;};


    /**
     * This function applies the local Stiffnessmatrices on all cells of  a fixed Depth.
     *
     *
     * @param input input data
     * @param output output data
     * @param depth
     */
    virtual void applyLocalStiffnessMatricesDepth(VectorSparseG& input, VectorSparseG& output, Depth& depth);

    void applyLocalStiffnessMatricesIndex(VectorSparseG& input, VectorSparseG& output, Depth& depth,IndexDimension& Index){
        //TODO
    };

    void applyLocalStiffnessMatricesFixedDepth(VectorSparseG& input, VectorSparseG& output, Depth& depth);
    void applyLocalStiffnessMatricesFixedDepth_onNode(VectorSparseG& input, VectorSparseG& output, Depth& D);

    double applyLocalStiffnessMatricesFixedDepthIndex_onNode(VectorSparseG& input, VectorSparseG& output, Depth& D, IndexDimension& Index);

    void applyLocalStiffnessMatricesOnIndices_onNode(VectorSparseG& input, VectorSparseG& output, Depth& D, std::vector<IndexDimension>& Indices);

    void receiveApplySend(int n);
    virtual void receiveApplySendOnActiveWorkers();

    void resetActiveWorkers();

    bool distribute=true;


    std::vector<int> active_worker;
    int getNumberProcesses(){return number_processes;};

    TypeMatrixVectorMultiplication getTypeMatrixVectorMultiplication(){return  typeMatrixVectorMultiplication;};

    int getNodeForDepth(Depth& T){return distributedDepthsHashtable.getNodeForDepth(T);};
    int numbercells=0;

    void releaseData(){
        pairedDepthsLocalStiffnessMatrices.clear();
    }

    int size(){ return  pairedDepthsLocalStiffnessMatrices.size();};



    double getTime(){
        double local_value = time;
        double global_max;
        MPI_Reduce(&local_value, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Bcast(&global_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


        if(abs(time-global_max)<1e-5){
            Tmax.Print();
            cout << numbercells_max << " " << time << "  " << global_max << endl;
        }

        local_value = time_min;
        double global_min;
        MPI_Reduce(&local_value, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Bcast(&global_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


        if(abs(time_min-global_min)<1e-5){

            cout << numbercells_min << " " << time_min << "  " << global_min << endl;
        }

        return global_max;
    }

    Depth getDepth(){

        return Tmax;};

    double time = 0.0;
    double time_min =1e+10;
    int numbercells_max=0;
    int numbercells_min=0;
    Depth Tmax;

    std::vector<std::pair<Depth,std::vector<LocalStiffnessMatrixFixedDepthSymmetric>>>*  getPairedDepthsLocalStiffnessMatrices(){
        return &pairedDepthsLocalStiffnessMatrices;
    }
protected:

    StencilTemplate& stencilClass;
    AdaptiveSparseGrid& grid;
    LocalStiffnessMatrices(AdaptiveSparseGrid& sg, StencilTemplate& stencilClass_,  int number_processes_):grid(sg), depthList(sg),stencilClass(stencilClass_), cellData(sg),input_node(sg),output_node(sg),number_processes(number_processes_){
/*        int rank;
        int size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);*/


    } ;

    int number_processes=1;
    void printEstimatedStorage();
    TypeMatrixVectorMultiplication typeMatrixVectorMultiplication = StoreLocalStiffnessMatrix;



    DepthList depthList;
    CellData cellData;
    DistributedDepthsHashtable distributedDepthsHashtable;

    VectorSparseG input_node,output_node;



    class LocalStiffnessMatrixFixedDepthSymmetric{
        friend LocalStiffnessMatrices;
    public:
        LocalStiffnessMatrixFixedDepthSymmetric( AdaptiveSparseGrid& grid_,StencilTemplate& stencilTemplate1):grid(grid_),stencilTemplate(stencilTemplate1){ };
        LocalStiffnessMatrixFixedDepthSymmetric(CellDimension& cell, AdaptiveSparseGrid& grid_, StencilTemplate& stencilTemplate_) :stencilTemplate(stencilTemplate_),grid(grid_){
            setValues(cell, stencilTemplate_);
            cellDimension=cell;

        }

        // Copy assignment operator
        LocalStiffnessMatrixFixedDepthSymmetric& operator=(const LocalStiffnessMatrixFixedDepthSymmetric& other) {
            if (this == &other) {
                return *this; // Handle self-assignment
            }
            for(int i=0; i <TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value; i++ ){
                occupied[i]=other.occupied[i];
                entries[i]=other.entries[i];
                storage_q[i]=other.storage_q[i];
                storage_p[i]=other.storage_p[i];
            }
            return *this;
        }


        // add operator
       void operator+=(const LocalStiffnessMatrixFixedDepthSymmetric& other) {

            for(int i=0; i <TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value; i++ ){
                if(occupied[i]) {
                    entries[i] += other.entries[i];

                    if (storage_q[i] != other.storage_q[i]) {
                        cout << " LocalStiffnessMatrixFixedDepthSymmetric : added wrong values " << storage_q[i] << "  "
                             << other.storage_q[i] << endl;
                        std::terminate();
                    }
                }
            }

        }

        bool operator==(const LocalStiffnessMatrixFixedDepthSymmetric& other) {

            bool retbool=true;
            for(int i=0; i <TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value; i++ ){
                if(occupied[i]) {
                    if(entries[i] != other.entries[i]){

                        cellDimension.PrintCoordinates();

                        retbool=false;
                        break;
                    }
                }
            }
            return retbool;

        }

        CellDimension* getCell(){return &cellDimension;};


        int array_size=TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value;   // Declaration of static member

        bool occupied[TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value];
        double entries[TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value];
        unsigned long storage_p[TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value];
        unsigned long storage_q[TriangularNumber<PowerOfTwo<DimensionSparseGrid>::value>::value];




        AdaptiveSparseGrid& grid;
        StencilTemplate& stencilTemplate;
        CellDimension cellDimension;



        void setValues(CellDimension& cell);

        void setValues(CellDimension& cell, StencilTemplate& stencil);

        void applyLocalMatrix(VectorSparseG& input, VectorSparseG& output);

        void applyLocalMatrixIndex(VectorSparseG& input, VectorSparseG& output,IndexDimension& Index);

        bool indexInMatrix(IndexDimension& Index);;



    };


/**
* @brief A vector containing pairs of Depth and arrays of doubles.
*
* This vector holds pairs, each consisting of a Depth object and a vector
* of LocalStiffnessMatrices for each cell.
*/
std::vector<std::pair<Depth,std::vector<LocalStiffnessMatrixFixedDepthSymmetric>>> pairedDepthsLocalStiffnessMatrices;
};















#endif //RUN_LOCALSTIFFNESSMATRICES_H
