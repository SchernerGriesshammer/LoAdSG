//
// Created by to35jepo on 12/7/22.
//

#ifndef SGRUN_MATRIXVECTORHOMOGENH_H
#define SGRUN_MATRIXVECTORHOMOGENH_H

#include "../extemp/vector.h"
#include "../tests/old_versions/MatrixVectorMultiplicationPrewavelets.h"
#include "../sgrid/multiDepthHashGrid.h"
#include "../BasisTransformations/BasisTransformations.h"
#include "../stencils/PoissonStencil.h"
#include "../stencils/Stencil.h"
#include "../iterator/depthIterator.h"
#include "../cells/CellStructure.h"
#include "../applications/norms.h"
#include "../mycomm.h"
#include "../localStiffnessMatrices/LocalStiffnessMatrices.h"



#include <omp.h>
#include <thread>
#include <sstream>
#include <sys/time.h>



/**
 *  This class provides a Matrix-Vector-Multiplication on Dirichlet-SparseGrids
 */
class MatrixVectorHomogen {
public:
    std::vector<Depth> depths;


    /**
     * @brief Constructor for MatrixVectorHomogen.
     * @param grid AdaptiveSparseGrid, make sure that the grid is completed, i.e. fullfills Definition of locally adaptive sparse grid
     * @param mgrid locally adaptive full subgrids
     */
     MatrixVectorHomogen(AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid) : gM(mgrid), g(grid), z(grid), u(grid),
                                                                                         nodal(mgrid),
                                                                                         u_new(grid),
                                                                                         u_old(grid), P(mgrid), Ax_neu(grid), nodal_test(mgrid), Ax_solution(grid),
                                                                                         depthList(grid), cellData(grid){};


    bool option_symmetry = false; /**< make use of symmetry of the bilinear form */
    bool option_mpi = false /**< Distribute local stiffness matrices with MPI. (No speedup yet!) */;


    /**
    * @brief Constructor for MatrixVectorHomogen.
    * @param grid AdaptiveSparseGrid, make sure that the grid is completed, i.e. fullfills Definition of locally adaptive sparse grid
    * @param mgrid locally adaptive full subgrids
     *@param symmetry make use of symmetry of the bilinear form
    */
    MatrixVectorHomogen(AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid, bool symmetry) : gM(mgrid), g(grid), z(grid), u(grid),
                                                                                         nodal(mgrid),
                                                                                         u_new(grid),
                                                                                         u_old(grid), P(mgrid), Ax_neu(grid), nodal_test(mgrid), Ax_solution(grid),
                                                                                         depthList(grid), option_symmetry(
                    symmetry),cellData(grid){};
    /**
    * @brief Constructor for MatrixVectorHomogen.
    * @param grid AdaptiveSparseGrid, make sure that the grid is completed, i.e. fullfills Definition of locally adaptive sparse grid
    * @param mgrid locally adaptive full subgrids
     *@param symmetry make use of symmetry of the bilinear form
     *@param mpi_on distribute local stiffness matrices with MPI
    */
    MatrixVectorHomogen(AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid, bool symmetry, bool mpi_on) : gM(mgrid), g(grid), z(grid), u(grid),
                                                                                                        nodal(mgrid),
                                                                                                        u_new(grid),
                                                                                                        u_old(grid), P(mgrid), Ax_neu(grid), nodal_test(mgrid), Ax_solution(grid),
                                                                                                        depthList(grid), option_symmetry(
                    symmetry),option_mpi(mpi_on),cellData(grid){};


    int stencil_count=0;/**< counter for stencil evaluations, just needed for testing */;

    /**
     * Calculates \f$ a(prew, \varphi_p) \quad \forall p \in D \f$ and stores it in Ax, with prew given in prewavelet basis.
     *
     *
     * @tparam StencilClass
     * @param prew
     * @param Ax
     * @param stencilClass
     */
    template<class StencilClass>
    void multiplication(VectorSparseG &prew, VectorSparseG &Ax, StencilClass &stencilClass);


    void multiplication_LocalStiffnessMatrix(VectorSparseG &prew, VectorSparseG& Ax, LocalStiffnessMatrices& localStiffnessMatrixAllDepths_){


        Ax = 0.0;
        nodal = 0.0;

        AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();
        int world_rank = 0;




        grid->WorkOnHangingNodes = true;




#pragma omp parallel
        {

            MultiDepthHashGrid MultiDepthHashGrid = *grid->getMultiDepthHashGrid();
            VectorSparseG u_nodal(u.getSparseGrid());
            int i = 0;

            for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {

                if (i % omp_get_num_threads() == omp_get_thread_num()) {

                    Depth T = *it;
                    u_nodal = 0.0;
                    calcNodal(u_nodal, prew, T);
                    nodal.setMultiLevelValuesOMP(u_nodal, T, MultiDepthHashGrid);
                }
                i++;
            }


        }



        grid->WorkOnHangingNodes = false;




        stencil_count = 0;

        for (CasesIterator iter2; iter2.goon(); ++iter2) {

            Ax_neu = 0.0;

            bool *restrictions = iter2.getcase();
            grid->WorkOnHangingNodes = true;


            Ax_neu = 0.0;
            gM = 0.0;


            stencil_count = 0;

            caseFunction_LocalStiffnessMatrix(prew, restrictions, localStiffnessMatrixAllDepths_);

            grid->WorkOnHangingNodes = false;
            Ax = Ax + Ax_neu;

        }




};

    template<class Problem>
    void multiplication_LocalStiffnessMatrix_mpi(VectorSparseG &prew, VectorSparseG& Ax, Problem& localStiffnessMatrixAllDepths_){

        localStiffnessMatrixAllDepths_.resetActiveWorkers();

#ifdef MY_MPI_ON
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        Ax = 0.0;
        nodal = 0.0;

        AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();





        grid->WorkOnHangingNodes = true;




#pragma omp parallel
        {

            MultiDepthHashGrid MultiDepthHashGrid = *grid->getMultiDepthHashGrid();
            VectorSparseG u_nodal(u.getSparseGrid());
            int i = 0;

            for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {

                if (i % omp_get_num_threads() == omp_get_thread_num()) {

                    Depth T = *it;
                    u_nodal = 0.0;
                    calcNodal(u_nodal, prew, T);
                    nodal.setMultiLevelValuesOMP(u_nodal, T, MultiDepthHashGrid);
                }
                i++;
            }


        }



        grid->WorkOnHangingNodes = false;



#ifdef MY_MPI_ON

        MPI_Barrier(MPI_COMM_WORLD);
        casefunction_mpi_onlyCases_LocalStiffnessMatrix(prew,Ax,localStiffnessMatrixAllDepths_);


         MPI_Barrier(MPI_COMM_WORLD);
#else

        stencil_count = 0;

        for (CasesIterator iter2; iter2.goon(); ++iter2) {

            Ax_neu = 0.0;

            bool *restrictions = iter2.getcase();
            grid->WorkOnHangingNodes = true;


            Ax_neu = 0.0;
            gM = 0.0;


            stencil_count = 0;

            caseFunction_LocalStiffnessMatrix(prew, restrictions, localStiffnessMatrixAllDepths_);

            grid->WorkOnHangingNodes = false;
            Ax = Ax + Ax_neu;

        }
#endif


#ifdef MY_MPI_ON
        MPI_Barrier(MPI_COMM_WORLD);
#endif

    };


    int mpi_per_worker=1;
    template<class StencilClass>
    void casefunction_mpi(VectorSparseG &prew, VectorSparseG &Ax, StencilClass &stencilClass);



    template<class StencilClass>
    void casefunction_mpi_onlyCases(VectorSparseG &prew, VectorSparseG &Ax, StencilClass &stencilClass);


    template<class Problem>
    void casefunction_mpi_onlyCases_LocalStiffnessMatrix(VectorSparseG &prew, VectorSparseG &Ax, Problem& localStiffnessMatrices);









    DepthList depthList;

private:

    MultiLevelVector gM;
    MultiLevelVector nodal;
    MultiLevelVector P;
    VectorSparseG Ax_neu;
    VectorSparseG g;
    VectorSparseG z;
    VectorSparseG u;
    VectorSparseG u_new;
    VectorSparseG u_old;
    MultiLevelVector nodal_test;
    VectorSparseG Ax_solution;

    CellData cellData;





    /**
     * Calculates nodal values of a function given of prewavelets (prew) on a fixed level (depth)
     * @param prew
     * @param depth
     */
    void calcNodal(VectorSparseG &prew, Depth &depth);

    /**
 * Calculates nodal values of a function given of prewavelets (prew) on a fixed level (depth)
 * @param prew
 * @param depth
 */
    void calcNodal(VectorSparseG &nodal_u,VectorSparseG &prew, Depth &depth);


    double time_sortdepth;
    double time_plus;

    /**
     * Calculates one of the 2^d cases.
     *
     *
     * @tparam StencilClass Stencil shall be given by a class
     * @param prew Function given in prewavelet basis
     * @param restrictions 2^d cases of prolongations (=0) and restrictions=1)
     * @param stencil
     */
    template<class StencilClass>
    void caseFunction(VectorSparseG &prew, bool *restrictions, StencilClass &stencil);



    /**
 * Calculates one of the 2^d cases.
 *
 *
 *
 * @param prew Function given in prewavelet basis
 * @param restrictions 2^d cases of prolongations (=0) and restrictions=1)
 * @param stencil
 */
    template<class Problem>
    void caseFunction_LocalStiffnessMatrix(VectorSparseG &prew, bool *restrictions, Problem& localStiffnessMatrixAllDepths);



    /**
     *
     * Apply a stencil on all grid points with depth less or equal to Depth  T
     *
     * @tparam StencilClass
     * @param input Nodal values
     * @param output a(input, v_p)
     * @param stencil
     * @param T
     */
    template<class StencilClass>
    void applyStencilGlobal(const VectorSparseG &input, VectorSparseG &output, StencilClass &stencil, Depth &T);

    /**
     *
     * Apply a stencil on all grid points with depth less or equal to Depth  T. Make use of Symmetry of the bilinearform by calculating local stiffness matrices on cells.
     *
     * @tparam StencilClass
     * @param input Nodal values
     * @param output a(input, v_p)
     * @param stencil
     * @param T
     */
    template<class StencilClass>
    void applyStencilGlobal_symmetry(VectorSparseG &input, VectorSparseG &output, StencilClass &stencil, Depth &T);

    /**
 *
 * Apply a stencil on all grid points with depth less or equal to Depth  T. Make use of Symmetry of the bilinearform by calculating local stiffness matrices on cells.
 * Distribute the local stiffness matrices with MPI. (No SpeedUp yet.)
 * @tparam StencilClass
 * @param input Nodal values
 * @param output a(input, v_p)
 * @param stencil
 * @param T
 */
    template<class StencilClass>
    void applyStencilGlobal_MPI(VectorSparseG &input, VectorSparseG &output, StencilClass &stencil, Depth &T);


    /**
*
* Apply a stencil on all grid points with depth less or equal to Depth  T. Make use of Symmetry of the bilinearform by calculating local stiffness matrices on cells.
* Distribute the local stiffness matrices with MPI. (No SpeedUp yet.)
* @tparam StencilClass
* @param input Nodal values
* @param output a(input, v_p)
* @param stencil
* @param T
*/
    template<class StencilClass>
    void applyStencilGlobal_MPI_2(VectorSparseG &input, VectorSparseG &output, StencilClass &stencil, Depth &T);


    /**
*
* Apply a stencil on all grid points with depth less or equal to Depth  T. Make use of Symmetry of the bilinearform by calculating local stiffness matrices on cells.
* Distribute the local stiffness matrices with MPI. (No SpeedUp yet.)
* @tparam StencilClass
* @param input Nodal values
* @param output a(input, v_p)
* @param stencil
* @param T
*/
    template<class StencilClass>
    void applyStencilGlobal_MPI_3(VectorSparseG &input, VectorSparseG &output, StencilClass &stencil, Depth &T);






    /**
     *
     * Apply a stencil on a given grid point (Index) w.r.t. Depth T. I.e. a(u,v_p) with u = nodal, p = Index.
     *
     * @tparam StencilClass
     * @param Index grid point
     * @param input nodal values
     * @param stencil
     * @param T
     * @return
     */
    template<class StencilClass>
    double applyStencilLocal(const IndexDimension &Index, const VectorSparseG &input, StencilClass &stencil, Depth &T);


/**
 *
 * Prolongation operator --> prolongation on primal space
 *
 * @param values
 * @param t depth
 * @param d direction
 */
    inline void prolongation1D_inplace(VectorSparseG &values, int &t, int &d) {



        AdaptiveSparseGrid_Base *grid = values.getSparseGrid();

        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it){
            Depth Tlocal = *it;
            if (Tlocal.at(d) == t) {
                SingleDepthHashGrid& depthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
                const auto& mapping = depthGrid._mapPosToGridPos;
                if(depthGrid.getNumberOfEntries()>0) {
                    #pragma omp parallel for
                    for (size_t i = 0; i < mapping.size(); i++) {

                        IndexDimension Index = depthGrid._map.getIndexOfTable(i);

                        unsigned long j = mapping[i];
                        double value = values.getValue(j);
                        if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d)))
                            value = value + 0.5 * (values.getValue(Index.nextRight(d)) + values.getValue(Index.nextLeft(d)));

                        values.setValue(j, value);


                    }
                }
            }

        }
    }

/**
 *
 * Restriction operator -> Restriction of dual space
 *
 *
 * @param values
 * @param t
 * @param d
 */
    inline void restriction1D_inverted_inplace(const VectorSparseG &values, int t, int d) {

        AdaptiveSparseGrid_Base *grid = values.getSparseGrid();

        //double value;
        double valueL;
        double valueR;

        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {
            Depth Tlocal = *it;

            // SingleDepthHashGrid& coarseDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
            if (Tlocal.at(d) == t + 1 && Tlocal > 0) {
                SingleDepthHashGrid &depthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
                const auto &mapping = depthGrid._mapPosToGridPos;
                if (depthGrid.getNumberOfEntries() > 0) {

                    for (size_t i = 0; i < mapping.size(); i++) {
                        IndexDimension Index = depthGrid._map.getIndexOfTable(i);

                        unsigned long j = mapping[i];
                        double value = 0.5 * values.getValue(j);

                        if (t >= 0) {
                            IndexDimension rightIndex = Index.nextRight(d);
                            if (!rightIndex.isAtRightBoundary(d)) {
                                // values.getValue(rightIndex,value,singleDepthHashGrid)
//#pragma omp critical
                                values.addToValue(rightIndex, value);
                                // value += values.getValue(rightIndex);
                            }

                            IndexDimension leftIndex = Index.nextLeft(d);
                            if (!leftIndex.isAtRightBoundary(d)) {
                                // values.getValue(rightIndex,value,singleDepthHashGrid)
//#pragma omp critical
                                values.addToValue(leftIndex, value);
                                // value += values.getValue(rightIndex);
                            }
                        }
                    }
                }
            }

        }
    };

    /**
     *
     * Convert dualspace to prewavlet basis i.e. a(u,nodal_p) -> a(u,prewavelet_p)
     *
     * @param Ax_nodal a(u,nodal_p)
     * @param Ax a(u,prewavlet_p)
     * @param T T(p) \leq T \forall p in D
     */
    inline void ConvertToPrewavelet2(VectorSparseG &Ax_nodal, VectorSparseG &Ax, Depth &T) {
        SingleDepthHashGrid& depthGrid = Ax_nodal.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
        const auto& mapping = depthGrid._mapPosToGridPos;

#pragma omp parallel for schedule(runtime)
        for (size_t i = 0; i < mapping.size(); i++)
        {
            if(Ax_nodal.getSparseGrid()->getActiveTable()[mapping[i]]) {
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                double coeff = 0.0;
                IndexDimension J;

                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                    double basis_coeff = 1.0;
                    J = I.nextFiveP(&mc, T, &basis_coeff);
                    if (mc.goon())
                        coeff = coeff + Ax_nodal.getValue(J) * basis_coeff;
                }

                Ax.setValue(mapping[i], coeff);
            }
        }
    }

    /**
     * This function is used for dynamic scheduling of the 2^d different cases of prolongations and restrictions.
     * The master sends a case (bool restrictions) to the worker. The worker then starts calculation the required partial sum of the bilinear form.
     *
     * @param num_workers
     * @param total_tasks
     */
    void master(int num_workers,int total_tasks);


    int MPI_task_index=0;
    int MPI_total_tasks=PowerOfTwo<DimensionSparseGrid>::value;
    int MPI_num_workers=0;


    /**
 * This function is used for dynamic scheduling of the 2^d different cases of prolongations and restrictions.
 * The master sends a case (bool restrictions) to the worker. The worker then starts calculation the required partial sum of the bilinear form.
 *
 * @param num_workers
 * @param total_tasks
 */
    void master_start(int num_workers);

    void master_start_onlyCases(int num_workers);

    template<class Problem>
    void master_start_onlyCases_LS(int num_workers,Problem &localStiffnessMatrices);




    /**
 * This function is used for dynamic scheduling of the 2^d different cases of prolongations and restrictions.
 * The master sends a case (bool restrictions) to the worker. The worker then starts calculation the required partial sum of the bilinear form.
 *
 */
    void master_distribute();


    void master_distribute_onlyCases();

    template<class Problem>
    void master_distribute_onlyCases_LS(Problem& localStiffnessMatrices);

    /**
* This function is used for dynamic scheduling of the 2^d different cases of prolongations and restrictions.
* The master sends a case (bool restrictions) to the worker. The worker then starts calculation the required partial sum of the bilinear form.
*
*/
    void master_end_distribute();

    template<class LS_Matrix>
    void master_end_distribute_LocalStiffnessMatrix(LS_Matrix& localStiffnessMatrices);


    void distributeAndApplyLocalStiffnessMatrix(LocalStiffnessMatrices &localStiffnessMatrices);

    /**
     *
     * worker function for dynamic mpi scheduling of the cases
     *
     * @tparam StencilClass
     * @param prew
     * @param Ax
     * @param stencilClass
     */
    template<class StencilClass>
    void worker(VectorSparseG &prew, VectorSparseG &Ax,StencilClass &stencilClass);

    template<class Problem>
    void worker_localStiffnessMatrices(VectorSparseG &prew, VectorSparseG &Ax,Problem &localStiffnessMatrices);



    /**
     * Converts an integer to a binary and stores the resulting binary as a boolean.
     *
     * @param num
     * @param binary
     */
    inline void intToBinary(int num, bool binary[]) {
        for (int i = DimensionSparseGrid - 1; i >= 0; --i) {
            binary[i] = num & 1;
            num >>= 1;
        }
    }

};



/////////  INLINE TEMPLATE FUNCTIONS//////////////////////////////////////////////////////////


template<class Problem>
void MatrixVectorHomogen::multiplication(VectorSparseG &prew, VectorSparseG &Ax, Problem &matrixProblem) {


if(matrixProblem.getTypeMatrixVectorMultiplication()==StencilOnTheFly) {




    Ax = 0.0;
    nodal = 0.0;

    AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();
    int world_rank = 0;


#ifdef MY_MPI_ON
    // Process *mpi = grid->mpi;
    int k = 0;
    int num_cases = POW(2,DimensionSparseGrid);
    // Get the rank of the process

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // TODO use mpi->getMyRank() instead. But is not correctly implemented
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // TODO could be done better with own defined COMM but maybe not needed

    MPI_Barrier(MPI_COMM_WORLD); // TODO do we need this barrier?

#endif


    grid->WorkOnHangingNodes = true;


#pragma omp parallel
    {

        MultiDepthHashGrid MultiDepthHashGrid = *grid->getMultiDepthHashGrid();
        VectorSparseG u_nodal(u.getSparseGrid());
        int i = 0;

        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {

            if (i % omp_get_num_threads() == omp_get_thread_num()) {

                Depth T = *it;
                u_nodal = 0.0;
                calcNodal(u_nodal, prew, T);
                nodal.setMultiLevelValuesOMP(u_nodal, T, MultiDepthHashGrid);
            }
            i++;
        }


    }


    grid->WorkOnHangingNodes = false;


#ifdef MY_MPI_ON

    casefunction_mpi_onlyCases(prew,Ax,matrixProblem);


#else


    stencil_count = 0;

    for (CasesIterator iter2; iter2.goon(); ++iter2) {

        Ax_neu = 0.0;

        bool *restrictions = iter2.getcase();
        grid->WorkOnHangingNodes = true;


        Ax_neu = 0.0;
        gM = 0.0;


        stencil_count = 0;

        caseFunction(prew, restrictions, matrixProblem);

        grid->WorkOnHangingNodes = false;
        Ax = Ax + Ax_neu;

    }


#endif
}else{

    multiplication_LocalStiffnessMatrix_mpi(prew, Ax,matrixProblem);
    MPI_Barrier(MPI_COMM_WORLD);
}

}


/**
 * @brief Distribute the cases with dynamic mpi_scheduling
 *
 * @tparam Problem
 * @param prew
 * @param Ax
 * @param stencilClass
 */

template<class Problem>
void MatrixVectorHomogen::casefunction_mpi(VectorSparseG &prew, VectorSparseG &Ax, Problem &stencilClass) {

#ifdef MY_MPI_ON
    int rank, num_procs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    MPI_task_index=0;

    if (rank < mpi_per_worker) {


        //distribute the first tasks to the workers
        if(rank ==0) {

            master_start(num_procs);


            //start working
            while(MPI_task_index<MPI_total_tasks) {


                for(int j=1; j < mpi_per_worker; j++){
                   MPI_Send(&MPI_task_index, 1, MPI_INT, j, TAG_TASK, MPI_COMM_WORLD);
                }


                prew.getSparseGrid()->WorkOnHangingNodes = true;
                Ax_neu = 0.0;
                gM = 0.0;
                bool restrictions[DimensionSparseGrid];
                intToBinary(MPI_task_index, restrictions);
                MPI_task_index++;

                caseFunction(prew, restrictions, stencilClass);
                prew.getSparseGrid()->WorkOnHangingNodes = false;
                Ax = Ax + Ax_neu;


                for(int j=1; j < mpi_per_worker; j++){
                        int old_index;

                        MPI_Recv(&old_index, 1, MPI_INT, j, TAG_TASK, MPI_COMM_WORLD,&status);
                }
            }



            for(int j=1; j < mpi_per_worker; j++){
                int old_index;
                MPI_Send(&MPI_task_index, 1, MPI_INT, j, TAG_TERMINATE, MPI_COMM_WORLD);
            }
        master_end_distribute();

        }else{

            while(1){

                MPI_Recv(&MPI_task_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


                if (TAG_TERMINATE == status.MPI_TAG) {
                     break;
                }



                prew.getSparseGrid()->WorkOnHangingNodes = true;
                Ax_neu = 0.0;
                gM = 0.0;
                bool restrictions[DimensionSparseGrid];
                intToBinary(MPI_task_index, restrictions);

                caseFunction(prew, restrictions, stencilClass);
                prew.getSparseGrid()->WorkOnHangingNodes = false;
                Ax = Ax + Ax_neu;

                MPI_Send(&MPI_task_index, 1, MPI_INT, 0, TAG_TASK, MPI_COMM_WORLD);

            }





        }





    } else {

        // Worker processes
        worker(prew, Ax, stencilClass);
    }




    if(rank%mpi_per_worker==0 && rank > 0 ){
        int tag=0;
        MPI_Send(Ax.getDatatableVector(), Ax.getSparseGrid()->getMaximalOccupiedSecondTable(),MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }

    if(rank==0){

        for(int rank_j=mpi_per_worker; rank_j< num_procs; rank_j+=mpi_per_worker){
            Ax_neu=0.0;
            MPI_Recv(Ax_neu.getDatatableVector(), Ax_neu.getSparseGrid()->getMaximalOccupiedSecondTable(),MPI_DOUBLE, rank_j, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            Ax = Ax+Ax_neu;
        }
    }

#endif



};


template<class Problem>
void MatrixVectorHomogen::casefunction_mpi_onlyCases(VectorSparseG &prew, VectorSparseG &Ax, Problem &stencilClass) {

#ifdef MY_MPI_ON
    int rank, num_procs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    MPI_task_index=0;
    MPI_num_workers=min(num_procs,MPI_total_tasks);





        //distribute the first tasks to the workers
        if(rank ==0) {


            master_start_onlyCases(num_procs);



            //start working
            while(MPI_task_index<MPI_total_tasks) {

            prew.getSparseGrid()->WorkOnHangingNodes = true;
            Ax_neu = 0.0;
            gM = 0.0;
            bool restrictions[DimensionSparseGrid];
            intToBinary(MPI_task_index, restrictions);
            MPI_task_index++;

            caseFunction(prew, restrictions, stencilClass);

            prew.getSparseGrid()->WorkOnHangingNodes = false;
            Ax = Ax + Ax_neu;

           }



        master_end_distribute();

        }else if(rank<MPI_num_workers){

            // Worker processes

             worker(prew, Ax, stencilClass);

        }


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, Ax.getDatatableVector(),Ax.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);



#endif



};


template<class StencilClass>
void MatrixVectorHomogen::caseFunction(VectorSparseG &prew, bool *restrictions, StencilClass &stencil) {



    gM= 0.0;
    AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();


    grid->WorkOnHangingNodes = true;







    // Sortiere Tiefen nach Restriktionen/Prolongationen
    depthList.sortDepths(restrictions);



    std::list<Depth>* sortedDepths = depthList.getSortierteTiefen();
    if(depthList.getSortierteTiefen()->size()==0) return;




    // Prolongiere die nodalen Werte
    P += nodal;



   for(int d=0; d < DimensionSparseGrid; d++) {
       if(restrictions[d]==0){
           for (Depth T: *sortedDepths) {

               int t = T.at(d) + 1;
               Depth T_fine = T;
               T_fine.set(t, d);


               if(depthList.isIncluded(T_fine)){


                   u_new = 0.0;
                   u_new.setMultiLevelValues2(P,T);
                   prolongation1D_inplace(u_new, t, d);


                   P.addMultiLevelValues(u_new, T_fine);

               }
           }
       }
   }



    for (Depth T: *sortedDepths){

         grid->WorkOnHangingNodes = true;


        u= 0.0;
        u.setMultiLevelValues2(P,T);

        z = 0.0;


        stencil.initialize(T);




        if(option_mpi){
           applyStencilGlobal_MPI_3(u,z,stencil,T);

       }else {
           if (option_symmetry){
           applyStencilGlobal_symmetry(u, z, stencil, T);

           }
           else applyStencilGlobal(u, z, stencil, T);
       }

#ifdef MY_MPI_ON

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank==0){
            master_distribute_onlyCases();
        }

#endif

       z.addMultiLevelValues2(gM,T);



       Depth Tcoarse(T);
       for (int d = 0; d < DimensionSparseGrid; d++) {
           if (restrictions[d]) {
               Tcoarse.set(T.at(d) - 1, d);
               if(depthList.isIncluded(Tcoarse)){
                   restriction1D_inverted_inplace(z,Tcoarse.at(d),d);
                   g = z;
                   z.addMultiLevelValues2(gM,Tcoarse);
                   gM.setMultiLevelValues2(g, Tcoarse);
               }
           }
       }



       grid->WorkOnHangingNodes=false;

       if(depthList.isIncluded(Tcoarse))   ConvertToPrewavelet2(z, Ax_neu, Tcoarse);



    }




};






template<class Problem>
void MatrixVectorHomogen::applyStencilGlobal(const VectorSparseG &input, VectorSparseG &output, Problem &matrixProblem, Depth &T) {


    AdaptiveSparseGrid_Base *grid = input.getSparseGrid();

    for(unsigned long i=0;i<grid->getMaximalOccupiedSecondTable();i++){
        IndexDimension I = grid->getIndexOfTable(i);
        Depth T_local(I);

        if(T_local<=T) {
            double val;
            val = applyStencilLocal(I, input, matrixProblem, T);
            output.setValue(i, val);
        }
    }


/*
    #pragma omp parallel
    {
        auto iter = GetNextSmallerDepthsIteratorInner(T);
        do {
            Depth Tlocal = *iter;

            const SingleDepthHashGrid &depthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
            const auto &mapping = depthGrid._mapPosToGridPos;
            const auto &map = depthGrid._map;
            const size_t end = mapping.size();


            if (depthGrid.getNumberOfEntries() > 0) {

                #pragma omp for schedule(runtime)
                for (size_t i = 0; i < mapping.size(); i++) {

                    double val;
                    IndexDimension I = map.getIndexOfTable(i);

                    val = applyStencilLocal(I, input, matrixProblem, T);



                    output.setValue(mapping[i], val);

                }
            }


        } while (iter.next());
    }
*/





};

template<class Problem>
void MatrixVectorHomogen::applyStencilGlobal_symmetry(VectorSparseG &input, VectorSparseG &output, Problem &matrixProblem, Depth &T) {
depths.push_back(T);





    AdaptiveSparseGrid_Base *grid = input.getSparseGrid();
    int rank = 0;

            Depth Tlocal = T;

            ++Tlocal;


            const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tlocal);
            const auto &map = depthGrid._map;
            const auto end = depthGrid.getNumberOfEntries();


#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < end; i++) {

        CellDimension cellDimension = map.getIndexOfTable(i);

        matrixProblem.applyStencilOnCell(cellDimension, input, output);
        //matrixProblem.applyStencilOnCell_MPI_OMP(cellDimension,input,output);


    }


};

template<class Problem>
void MatrixVectorHomogen::applyStencilGlobal_MPI(VectorSparseG &input, VectorSparseG &output, Problem &matrixProblem, Depth &T) {

    Depth Tlocal = T;
    ++Tlocal;
#ifdef MY_MPI_ON
    int rank; int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);




    int case_size = mpi_per_worker;
    int case_rank = rank%mpi_per_worker;

#endif


    const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tlocal);
    const auto &map = depthGrid._map;
    const auto end_grid = depthGrid.getNumberOfEntries();

        for (size_t i = 0; i < end_grid; i++) {


#ifdef  MY_MPI_ON
                if (i % case_size == case_rank) {

#endif

                    CellDimension cellDimension = map.getIndexOfTable(i);
                    matrixProblem.applyStencilOnCell_MPI_OMP(cellDimension, input, output);
#ifdef  MY_MPI_ON
            }
#endif

        }



#ifdef MY_MPI_ON



        //this should be done with MPI communicators
    if(case_rank==0){
            //case_rank == 0 means master

            VectorSparseG output_sum(output.getSparseGrid());


            for(int j=1;j<case_size;j++){

                MPI_Status status;


                MPI_Recv(output_sum.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, rank+j, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                output = output+output_sum;
             
             }

             for(int j=1;j<case_size;j++){


                 int tag=0;
                 MPI_Send(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE,rank+j, tag, MPI_COMM_WORLD);

             }




    }else{



            int sendTo = rank-case_rank;

            int tag=0;
            MPI_Status status;

            MPI_Send(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE,sendTo, tag, MPI_COMM_WORLD);

            MPI_Recv(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, sendTo, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    }




#endif




};




template<class Problem>
void MatrixVectorHomogen::applyStencilGlobal_MPI_2(VectorSparseG &input, VectorSparseG &output, Problem &matrixProblem, Depth &T) {

    Depth Tlocal = T;
    ++Tlocal;
#ifdef MY_MPI_ON
    int rank; int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);




    int case_size = mpi_per_worker;
    if(mpi_per_worker!=2){
        cout << "only 2 mpi workers are allowed " << endl;
        exit(1);
    }
    int case_rank = rank%mpi_per_worker;
#else
    cout << "THIS CAN ONLY USED WITH MPI! " << endl;


#endif


    const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tlocal);
    const auto &map = depthGrid._map;
    const auto end_grid = depthGrid.getNumberOfEntries();

    if(end_grid%2==0){
#ifdef  MY_MPI_ON
        if (case_rank==0) {
            #pragma omp parallel for schedule(runtime)
            for (size_t i = 0; i < end_grid/2; i++) {
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    matrixProblem.applyStencilOnCell(cellDimension, input, output);
            }
        }
        if (case_rank==1){
            #pragma omp parallel for schedule(runtime)
            for(size_t i=end_grid/2; i < end_grid; i++){
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    matrixProblem.applyStencilOnCell(cellDimension, input, output);
            }
        }
#endif
    }else{
#ifdef  MY_MPI_ON
        if (case_rank==0) {
            #pragma omp parallel for schedule(runtime)
            for (size_t i = 0; i < (end_grid-1)/2; i++) {
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    matrixProblem.applyStencilOnCell(cellDimension, input, output);
            }
        }
        if (case_rank==1){
            #pragma omp parallel for schedule(runtime)
            for(size_t i=(end_grid-1)/2; i < end_grid; i++){
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    matrixProblem.applyStencilOnCell(cellDimension, input, output);
            }
        }
#endif


    }




#ifdef MY_MPI_ON




    if(case_rank==0){


            VectorSparseG output_sum(output.getSparseGrid());
            MPI_Status status;
            MPI_Recv(output_sum.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            output = output+output_sum;
            int tag=0;
            MPI_Send(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE,rank+1, tag, MPI_COMM_WORLD);

    }else{
            int sendTo = rank-1;
            int tag=0;
            MPI_Status status;
            MPI_Send(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE,sendTo, tag, MPI_COMM_WORLD);
            MPI_Recv(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, sendTo, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    }




#endif




};




template<class Problem>
void MatrixVectorHomogen::applyStencilGlobal_MPI_3(VectorSparseG &input, VectorSparseG &output, Problem &matrixProblem, Depth &T) {

    Depth Tlocal = T;
    ++Tlocal;
#ifdef MY_MPI_ON
    int rank; int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);




    int case_size = mpi_per_worker;

    int case_rank = rank%mpi_per_worker;
#else
    cout << "THIS CAN ONLY USED WITH MPI! " << endl;


#endif


    const SingleDepthHashCellStructure &depthGrid = cellData.getGridForDepth(Tlocal);
    const auto &map = depthGrid._map;
    const auto end_grid = depthGrid.getNumberOfEntries();

    if(end_grid%mpi_per_worker==0) {
#ifdef  MY_MPI_ON

        for(int j=0; j<mpi_per_worker; j++){
            if (case_rank==j){
#pragma omp parallel for schedule(runtime)
                for (size_t i = j*(end_grid/mpi_per_worker); i < end_grid/mpi_per_worker; i++) {
                        CellDimension cellDimension = map.getIndexOfTable(i);
                        matrixProblem.applyStencilOnCell(cellDimension, input, output);
                }
            }
        }

#endif
    }else{
#ifdef  MY_MPI_ON
        int modulo = end_grid%mpi_per_worker;
        int end_grid_new = end_grid-modulo;
        for(int j=0; j<mpi_per_worker; j++){
            if (case_rank==j){
            #pragma omp parallel for schedule(runtime)
                for (size_t i = j*(end_grid_new/mpi_per_worker); i < end_grid_new/mpi_per_worker; i++) {
                    CellDimension cellDimension = map.getIndexOfTable(i);
                    matrixProblem.applyStencilOnCell(cellDimension, input, output);
                }
            }
        }

        for(int j=0; j<modulo; j++){
            if (case_rank==j){
                CellDimension cellDimension = map.getIndexOfTable(j+end_grid_new);
                matrixProblem.applyStencilOnCell(cellDimension, input, output);
            }
        }



#endif


    }




#ifdef MY_MPI_ON




    //this should be done with MPI communicators
    if(case_rank==0){
            //case_rank == 0 means master

            VectorSparseG output_sum(output.getSparseGrid());


            for(int j=1;j<case_size;j++){

                MPI_Status status;


                MPI_Recv(output_sum.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, rank+j, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                output = output+output_sum;

             }

             for(int j=1;j<case_size;j++){


                 int tag=0;
                 MPI_Send(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE,rank+j, tag, MPI_COMM_WORLD);

             }




    }else{



            int sendTo = rank-case_rank;

            int tag=0;
            MPI_Status status;

            MPI_Send(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE,sendTo, tag, MPI_COMM_WORLD);

            MPI_Recv(output.getDatatableVector(), output.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, sendTo, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    }





#endif




};




template<class Problem>
double MatrixVectorHomogen::applyStencilLocal(const IndexDimension &Index, const VectorSparseG &input, Problem &matrixProblem, Depth &T) {
    double val = 0.0;
    IndexDimension NextIndex;



    for (MultiDimCompass mc; mc.goon(); ++mc) {
        NextIndex = Index.nextThree2(&mc, T.returnTiefen());
        unsigned long k;


        if(mc.goon()&&input.getSparseGrid()->occupied(k,NextIndex)) {
            double value = 0.0;



            value = matrixProblem.returnValue(const_cast<IndexDimension &>(Index), mc);





            val = val + value * input.getValue(k);





        }


    }
    return val;

};



template<class Problem>
void MatrixVectorHomogen::worker(VectorSparseG &prew, VectorSparseG &Ax,Problem &matrixProblem) {
#ifdef MY_MPI_ON
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int task_index;
    MPI_Status status;

        while (1) {

            MPI_Recv(&task_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);


            if (TAG_TERMINATE == status.MPI_TAG) {

                break;  // Terminate worker
            }

            prew.getSparseGrid()->WorkOnHangingNodes = true;

            Ax_neu = 0.0;
            gM = 0.0;
            bool restrictions[DimensionSparseGrid];
            intToBinary(task_index, restrictions);


            caseFunction(prew, restrictions, matrixProblem);

            prew.getSparseGrid()->WorkOnHangingNodes = false;
            Ax = Ax + Ax_neu;


            MPI_Send(&task_index, 1, MPI_INT, 0, TAG_TASK, MPI_COMM_WORLD);


        }
#endif

}






/////////////////// TEMPLATED FUNCTIONS
template <class Problem>
void MatrixVectorHomogen::caseFunction_LocalStiffnessMatrix(VectorSparseG &prew, bool *restrictions, Problem& localStiffnessMatrixAllDepths) {
    int rank;



    gM= 0.0;
    AdaptiveSparseGrid_Base *grid = prew.getSparseGrid();


    grid->WorkOnHangingNodes = true;







    // Sortiere Tiefen nach Restriktionen/Prolongationen
    depthList.sortDepths(restrictions);



    std::list<Depth>* sortedDepths = depthList.getSortierteTiefen();
    if(depthList.getSortierteTiefen()->size()==0) return;




    // Prolongiere die nodalen Werte
    P += nodal;





    for(int d=0; d < DimensionSparseGrid; d++) {
        if(restrictions[d]==0){
            for (Depth T: *sortedDepths) {

                int t = T.at(d) + 1;
                Depth T_fine = T;
                T_fine.set(t, d);

                if(depthList.isIncluded(T_fine)){

                    u_new = 0.0;
                    u_new.setMultiLevelValues2(P,T);
                    prolongation1D_inplace(u_new, t, d);


                    P.addMultiLevelValues(u_new, T_fine);

                }
            }
        }
    }

    for (Depth T: *sortedDepths) {

        grid->WorkOnHangingNodes = true;


        u= 0.0;
        u.setMultiLevelValues2(P,T);

        z = 0.0;




        localStiffnessMatrixAllDepths.applyLocalStiffnessMatricesDepth(u,z,T);


#ifdef MY_MPI_ON

        int rank;
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                if(rank==0){


                    master_distribute_onlyCases_LS(localStiffnessMatrixAllDepths);
                }

#endif




        z.addMultiLevelValues2(gM,T);



        Depth Tcoarse(T);
        for (int d = 0; d < DimensionSparseGrid; d++) {
            if (restrictions[d]) {
                Tcoarse.set(T.at(d) - 1, d);
                if(depthList.isIncluded(Tcoarse)){
                    restriction1D_inverted_inplace(z,Tcoarse.at(d),d);
                    g = z;
                    z.addMultiLevelValues2(gM,Tcoarse);
                    gM.setMultiLevelValues2(g, Tcoarse);
                }
            }
        }



        grid->WorkOnHangingNodes=false;


        if(depthList.isIncluded(Tcoarse))   ConvertToPrewavelet2(z, Ax_neu, Tcoarse);



    }




}


template <class Problem>
void MatrixVectorHomogen::master_start_onlyCases_LS(int num_workers,Problem &localStiffnessMatrices) {
#ifdef MY_MPI_ON




    MPI_Status status;
    MPI_task_index=0;

    int num_tasks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);





    for (int worker_rank = 1; (worker_rank < num_workers-localStiffnessMatrices.getNumberProcesses()) && (MPI_task_index < MPI_total_tasks); worker_rank++) {

        MPI_Send(&MPI_task_index, 1, MPI_INT, worker_rank, TAG_TASK, MPI_COMM_WORLD);


        MPI_task_index++;
    }



#endif

}


template <class Problem>
void MatrixVectorHomogen::worker_localStiffnessMatrices(VectorSparseG &prew, VectorSparseG &Ax,
                                                        Problem &localStiffnessMatrices) {
#ifdef MY_MPI_ON

    int rank;
    int num_tasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);



    int task_index;
    MPI_Status status;
    int flag=0;
    if(rank<num_tasks-localStiffnessMatrices.getNumberProcesses()){

        while (1){

                   flag =0;
                   MPI_Iprobe(0,MPI_ANY_TAG,MPI_COMM_WORLD,&flag, &status);
                   if(flag){

                   MPI_Recv(&task_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                    if (TAG_TERMINATE == status.MPI_TAG) {
                        break;  // Terminate worker
                    }



                    prew.getSparseGrid()->WorkOnHangingNodes = true;

                    Ax_neu = 0.0;
                    gM = 0.0;
                    bool restrictions[DimensionSparseGrid];
                    intToBinary(task_index, restrictions);


                    caseFunction_LocalStiffnessMatrix(prew, restrictions, localStiffnessMatrices);

                    prew.getSparseGrid()->WorkOnHangingNodes = false;
                    Ax = Ax + Ax_neu;


                    MPI_Send(&task_index, 1, MPI_INT, 0, TAG_TASK, MPI_COMM_WORLD);
                    }

        }


        for (int n = num_tasks-localStiffnessMatrices.getNumberProcesses(); n < num_tasks; n++) {

            MPI_Send(&MPI_task_index, 1, MPI_INT, n, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD);
        }
    }else {


        while(localStiffnessMatrices.active_worker.size()>0){

            localStiffnessMatrices.receiveApplySendOnActiveWorkers();
        }
    }


#endif

}



template <class Problem>
void MatrixVectorHomogen::casefunction_mpi_onlyCases_LocalStiffnessMatrix(VectorSparseG &prew, VectorSparseG &Ax, Problem &localStiffnessMatrices) {

#ifdef MY_MPI_ON
    int rank, num_procs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    MPI_task_index=0;
    MPI_num_workers=num_procs;
    Ax = 0.0;





        //distribute the first tasks to the workers
        if(rank ==0){



            master_start_onlyCases_LS(num_procs,localStiffnessMatrices);



            //start working

            while(MPI_task_index<MPI_total_tasks) {


            prew.getSparseGrid()->WorkOnHangingNodes = true;
            Ax_neu = 0.0;
            gM = 0.0;
            bool restrictions[DimensionSparseGrid];
            intToBinary(MPI_task_index, restrictions);
            MPI_task_index++;

            caseFunction_LocalStiffnessMatrix(prew, restrictions, localStiffnessMatrices);

            prew.getSparseGrid()->WorkOnHangingNodes = false;
            Ax = Ax + Ax_neu;

            }


            master_end_distribute_LocalStiffnessMatrix(localStiffnessMatrices);

        }else{

            // Worker processes


            worker_localStiffnessMatrices(prew, Ax, localStiffnessMatrices);





        }


        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Allreduce(MPI_IN_PLACE, Ax.getDatatableVector(),Ax.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
/*
        int n=num_procs-localStiffnessMatrices.getNumberProcesses();
        int ranks_mv[n];
        for(int i=0;i<n;i++)ranks_mv[i]=i;

        MPI_Group world_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);
        MPI_Group group1;
        MPI_Group_incl(world_group, n, ranks_mv, &group1);
        MPI_Comm comm1;
        MPI_Comm_create(MPI_COMM_WORLD, group1, &comm1);


    if (comm1 != MPI_COMM_NULL)MPI_Allreduce(MPI_IN_PLACE, Ax.getDatatableVector(),Ax.getSparseGrid()->getMaximalOccupiedSecondTable(), MPI_DOUBLE, MPI_SUM,comm1);

    MPI_Group_free(&world_group);
    MPI_Group_free(&group1);

    if (comm1 != MPI_COMM_NULL) {
        MPI_Comm_free(&comm1);
    }
*/
    MPI_Barrier(MPI_COMM_WORLD);
#endif



};


template <class Problem>
void MatrixVectorHomogen::master_distribute_onlyCases_LS(Problem& localStiffnessMatrices) {
#ifdef MY_MPI_ON

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    MPI_Status status;
    int flag=0;
    int old_index;





    for (int worker =1; worker < MPI_num_workers-localStiffnessMatrices.getNumberProcesses() && MPI_task_index < MPI_total_tasks; worker++) {

            MPI_Iprobe(worker, TAG_TASK, MPI_COMM_WORLD, &flag, &status);

            if (flag){

                    MPI_Recv(&old_index, 1, MPI_INT, worker, TAG_TASK, MPI_COMM_WORLD, &status);

                    MPI_Send(&MPI_task_index, 1, MPI_INT, worker, TAG_TASK, MPI_COMM_WORLD);
                    MPI_task_index++;
            }


    }






#endif

};




template <class LS_Matrix>
void MatrixVectorHomogen::master_end_distribute_LocalStiffnessMatrix(LS_Matrix &localStiffnessMatrices){

#ifdef MY_MPI_ON



    MPI_Status status;
    int flag=0;
    int old_index;



        for (int n =  MPI_num_workers-localStiffnessMatrices.getNumberProcesses(); n <  MPI_num_workers; n++) {

            MPI_Send(&MPI_task_index, 1, MPI_INT, n, TAG_LOCALSTIFFNESS_END, MPI_COMM_WORLD);

        }


        // Send termination signal to workers
        for (int worker_rank = 1; worker_rank < MPI_num_workers-localStiffnessMatrices.getNumberProcesses(); worker_rank++) {

                   //receive message, that worker has completed all of its tasks
                   MPI_Recv(&old_index, 1, MPI_INT, worker_rank, TAG_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                   //send termination signal to worker

                   MPI_Send(&MPI_task_index, 1, MPI_INT, worker_rank, TAG_TERMINATE, MPI_COMM_WORLD);


        }






#endif

}





#endif //SGRUN_MATRIXVECTORHOMOGENH_H

















/////////////// OTHER PROLONGATIONS AND RESTRICTIONS: MIGHT BE USED LATER

/*    inline void prolongation1D(VectorSparseG &coarse, VectorSparseG &fine, int &t, int &d) {
        AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();
        unsigned long maxocc = grid->getMaximalOccupiedSecondTable();

        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {
            Depth Tlocal = iter.getDepth();
            if (Tlocal.at(d) == t) {
                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                do {
                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    double value = fine.getValue(i);
                    if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d)))
                        value = value + 0.5 * (fine.getValue(Index.nextRight(d)) + fine.getValue(Index.nextLeft(d)));
                    //coarse = value | Index;
                    coarse.setValue(i, value);

                } while (inneriter.next());
            }

        } while (iter.next());

*//*

        for (unsigned long i = 0; i < maxocc; i++) {
            IndexDimension Index = grid->getIndexOfTable(i);

            if (Index.getDepth(d) == t) {
                double value = fine.getValue(i);
*//*
*//*
                    if (Index.isAtLeftBoundary(d))
                        value = fine.getValue(Index);

                    if (Index.isAtRightBoundary(d))
                        value = fine.getValue(Index);
*//**//*


                if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d)))
                    value = value + 0.5 *( fine.getValue(Index.nextRight(d))+fine.getValue(Index.nextLeft(d)));
                //coarse.setValue(Index, value);
                coarse = value | Index;


            }


        }
*//*


    }*/


/*    void prolongation1D_inplace_inverted(VectorSparseG &vec,const int &t,const int &d) {
        AdaptiveSparseGrid_Base *grid = vec.getSparseGrid();
        unsigned long maxocc = grid->getMaximalOccupiedSecondTable();


        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {

            Depth Tlocal = iter.getDepth();
            //         Tlocal.Print();
            // cout << "Nodes: " <<endl;

            if (Tlocal.at(d) < t) {
                Depth searchDepth = Tlocal;
                searchDepth.set(t,d);
                //SingleDepthHashGrid& singleDepthHashGrid = grid->getMultiDepthHashGrid()->getGridForDepth(searchDepth);
                // Tlocal.Print();
                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();

                // int i = 0;
                do {
                    // i++;
                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    double value = vec.getValue(i) * 0.5;
                    IndexDimension rightIndex = Index.nextRight(d,t);
                    if (!rightIndex.isAtRightBoundary(d)){
                        // vec.getValue(rightIndex,value,singleDepthHashGrid)
                        vec.addToValue(rightIndex,value);
                        // value += vec.getValue(rightIndex);
                    }
                    IndexDimension leftIndex = Index.nextLeft(d,t);
                    if (!leftIndex.isAtLeftBoundary(d)){
                        vec.addToValue(leftIndex,value);
                        // value += vec.getValue(leftIndex);
                    }
                } while (inneriter.next());
                // cout << i << ", " << (1<<Tlocal.LoneNorm()-3)<< endl;
            }

            // cout << "iter" << endl;
        } while (iter.next());
        // cout << "end" << endl;


    }*/

/*    void prolongation1D_inplace_singleHash(VectorSparseG &vec,const int &t,const int &d) {
        AdaptiveSparseGrid_Base *grid = vec.getSparseGrid();
        unsigned long maxocc = grid->getMaximalOccupiedSecondTable();



        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {

            Depth Tlocal = iter.getDepth();
            //         Tlocal.Print();
            // cout << "Nodes: " <<endl;

            if (Tlocal.at(d) == t) {


                vector<SingleDepthHashGrid*> fineDepthGrids = grid->getMultiDepthHashGrid()->getGridsForDepthInDirection(Tlocal,d);
                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                do {
                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    double value = vec.getValue(i);
                    if ((!Index.isAtRightBoundary(d)) && (!Index.isAtLeftBoundary(d))){
                        auto right = Index.nextRight(d);
                        auto left = Index.nextLeft(d);
                        // auto* gridRight = fineDepthGrids[right.getDepth(d)];
                        // auto* gridLeft = fineDepthGrids[left.getDepth(d)];
                        // auto valRight = vec.getValue(right,gridRight);
                        // auto valLeft = vec.getValue(left,gridLeft);
                        value = value + 0.5 * (vec.getValue(right,fineDepthGrids[right.getDepth(d)]) + vec.getValue(left,fineDepthGrids[left.getDepth(d)]));
                        // value = value + 0.5 * (vec.getValue(Index.nextRight(d)) + vec.getValue(Index.nextLeft(d)));
                        // auto p= vec.getValue(right,fineDepthGrids[right.getDepth(d)]);
                    }
                    //coarse = value | Index;
                    vec.setValue(i, value);
                    // cout << Depth(Index.nextRight(d)).at(d)<<", \t" <<Depth(Index.nextLeft(d)).at(d) << endl;

                    // cout << "Left: " ;
                    //     Depth(Index.nextLeft(d)).Print();
                    // cout << "Right: " ;
                    //     Depth(Index.nextRight(d)).Print();
                } while (inneriter.next());


            }

            // cout << "iter" << endl;
        } while (iter.next());
        // cout << "end" << endl;


    }*/

/*    inline void prolongation(VectorSparseG &coarse, VectorSparseG &fine, Depth &Tcoarse, Depth &Tfine) {
        fine = coarse;
        int tcoarse;
        int tfine;

        for (int d = 0; d < DimensionSparseGrid; d++) {

            tcoarse = Tcoarse.at(d) + 1;
            tfine = Tfine.at(d);
            for (int t = tcoarse; t <= tfine; ++t) {
                prolongation1D(fine, fine, t, d);
            }
        }
    }*/



/*    inline void restriction1D_inverted_inplace( const VectorSparseG &fine,const VectorSparseG &coarse,int t, int d) {


        double value=0.0;

        AdaptiveSparseGrid_Base* grid = fine.getSparseGrid();

        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {
            Depth Tlocal = iter.getDepth();
            // SingleDepthHashGrid& coarseDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
            if(Tlocal.at(d)<=t){

                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                // Depth fineDepth = Depth(inneriter.getPoint());
                // SingleDepthHashGrid& fineDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Depth(inneriter.getPoint().nextRight(d,t+1)));
                do {

                    unsigned long i = inneriter.geti();
                    IndexDimension I = inneriter.getPoint();
                    double value = fine.getValue(i);
                    unsigned long k;

                    coarse.setValue(i,value);

                } while (inneriter.next());

            }

            if (Tlocal.at(d) == t + 1 && Tlocal>0) {

                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                // Depth fineDepth = Depth(inneriter.getPoint());
                SingleDepthHashGrid& fineDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Depth(inneriter.getPoint().nextRight(d,t+1)));
                do {

                    IndexDimension Index = inneriter.getPoint();
                    unsigned long i = inneriter.geti();
                    value = 0.5 * fine.getValue(i);

                        if (t >= 0) {

                            IndexDimension rightIndex = Index.nextRight(d);
                            if (!Index.isAtRightBoundary(d)) {
                                // vec.getValue(rightIndex,value,singleDepthHashGrid)
                                coarse.addToValue(rightIndex, value);
                                // value += vec.getValue(rightIndex);
                            }


                            IndexDimension leftIndex = Index.nextLeft(d);
                            if (!Index.isAtLeftBoundary(d)) {
                                // vec.getValue(rightIndex,value,singleDepthHashGrid)
                                coarse.addToValue(leftIndex, value);
                                // value += vec.getValue(rightIndex);
                            }
                        }

                } while (inneriter.next());
            }
        } while (iter.next());

    };*/

/*   inline void restriction1D_inverted2( const VectorSparseG &fine, VectorSparseG &coarse, int t, int d) {

       AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();

       double value;
       double valueL;
       double valueR;

       ListOfDepthOrderedSubgrids::iterator iter(list);
       iter.gotoBegin();
       do {
           Depth Tlocal = iter.getDepth();
           // SingleDepthHashGrid& coarseDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
           if (Tlocal.at(d) == t + 1 && Tlocal > 0) {
               SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
               inneriter.gotobegin();
               // Depth fineDepth = Depth(inneriter.getPoint());
               // SingleDepthHashGrid& fineDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Depth(inneriter.getPoint().nextRight(d,t+1)));
               do {

                   IndexDimension Index = inneriter.getPoint();
                   unsigned long i = inneriter.geti();
                   value = 0.5 * fine.getValue(i);

                   if (t >= 0) {
                       IndexDimension rightIndex = Index.nextRight(d);
                       if (!rightIndex.isAtRightBoundary(d)){
                           // vec.getValue(rightIndex,value,singleDepthHashGrid)
                           coarse.addToValue(rightIndex,value);
                           // value += vec.getValue(rightIndex);
                       }

                       IndexDimension leftIndex = Index.nextLeft(d);
                       if (!leftIndex.isAtRightBoundary(d)){
                           // vec.getValue(rightIndex,value,singleDepthHashGrid)
                           coarse.addToValue(leftIndex,value);
                           // value += vec.getValue(rightIndex);
                       }
                   }
               } while (inneriter.next());
           }
       } while (iter.next());

   };*/



/*    inline void restriction1D(const VectorSparseG &fine, VectorSparseG &coarse, int t, int d) {

        AdaptiveSparseGrid_Base *grid = fine.getSparseGrid();

        double value;
        double valueL;
        double valueR;

        ListOfDepthOrderedSubgrids::iterator iter(list);
        iter.gotoBegin();
        do {
            Depth Tlocal = iter.getDepth();
            // SingleDepthHashGrid& coarseDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Tlocal);
            if (Tlocal.at(d) <= t && Tlocal > 0) {
                SubgridFixedDepth::iterator inneriter(*iter.getSubgrid());
                inneriter.gotobegin();
                // Depth fineDepth = Depth(inneriter.getPoint());
                SingleDepthHashGrid& fineDepthGrid = grid->getMultiDepthHashGrid()->getGridForDepth(Depth(inneriter.getPoint().nextRight(d,t+1)));
                do {

                    IndexDimension Index = inneriter.getPoint();



                    unsigned long i = inneriter.geti();
                    value = fine.getValue(i);
                    valueL = 0.0;
                    valueR = 0.0;
                    if (t >= 0) {
                        if (!(Index.isAtRightBoundary(d))) valueL = 0.5 * fine.getValue(Index.nextRight(d, t + 1));
                        if (!(Index.isAtLeftBoundary(d))) valueR = 0.5 * fine.getValue(Index.nextLeft(d, t + 1));
                    }


                    value = value + valueL + valueR;
                    //coarse.setValue(Index, value);
                    coarse.setValue(i, value);



                } while (inneriter.next());
            }
        } while (iter.next());

    };*/


/*
void prolongation_inplace(VectorSparseG &vec, const Depth &Tcoarse, const Depth &Tfine) {
    int tcoarse;
    int tfine;



    for (int d = 0; d < DimensionSparseGrid; d++) {

        tcoarse = Tcoarse.at(d) + 1;
        tfine = Tfine.at(d);
        for (int t = tcoarse; t <= tfine; ++t) {
            // if(t = tfine)
            //     prolongation1D_inplace_inverted(vec, t, d);
            // else
            prolongation1D_inplace(vec, t, d);
        }
    }
}*/
