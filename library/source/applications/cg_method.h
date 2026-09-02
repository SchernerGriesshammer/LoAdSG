//
// Created by scherner on 26.07.21.
//

#ifndef GRUN_CG_METHOD_H
#define GRUN_CG_METHOD_H


#include <sys/time.h>
#include <unordered_set>
#include "../MatrixVectorMultiplication/MatrixVectorHomogen.h"
#include "../MatrixVectorMultiplication/MatrixVectorInhomogen.h"

template<class StencilClass>
class Preconditioning{
public:
    Preconditioning(VectorSparseG& z_, StencilClass stencilClass, bool constant = false) : z(z_), boundary(false)
    {
        if (!constant){
            if(stencilClass.getTypeMatrixVectorMultiplication()==StencilOnTheFly) precond3(stencilClass);
            else{
                precondLocalStiffnessMatrix(stencilClass);
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }else
        {
            int level = z.getSparseGrid()->getLevel();
            double val = pow(2.0,level);
            for (unsigned long k = 0; k < z.getSparseGrid()->getMaximalOccupiedSecondTable(); k++)
            {
                if(z.getSparseGrid()->getActiveTable()[k])
                {
                    z.setValue(k,val);
                }
            }

        }
    }



    void precondMV(StencilClass stencilClass){
        int world_rank;
        int world_size;
#ifdef MY_MPI_ON

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif



        VectorSparseG v(*z.getSparseGrid());
        VectorSparseG Av(*z.getSparseGrid());
        MultiLevelAdaptiveSparseGrid mgrid(z.getSparseGrid());
        MatrixVectorHomogen matrixVectorHomogen(*z.getSparseGrid(), mgrid, 1, 0);

        for (size_t i = 0; i < z.getSparseGrid()->getMaximalOccupiedSecondTable(); i++) {
            if (z.getSparseGrid()->getActiveTable()[i]) {

                IndexDimension Index = z.getSparseGrid()->getIndexOfTable(i);


                v.setValue(Index,1);
                matrixVectorHomogen.multiplication(v,Av,stencilClass);




                v = 0.0;
                double coeff=Av.getValue(i);
                if(world_rank==0) {
                    Index.PrintCoord();
                    cout << " " << coeff << endl;
                }
                z.setValue(i, coeff);

            }


        }



    }

    void precond1(StencilClass stencilClass){
#ifdef MY_MPI_ON
        int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif



        VectorSparseG v(*z.getSparseGrid());

        for (size_t i = 0; i < z.getSparseGrid()->getMaximalOccupiedSecondTable(); i++) {
            if (z.getSparseGrid()->getActiveTable()[i]) {
#ifdef MY_MPI_ON
                if (world_rank == i % world_size ) {

#endif
                IndexDimension Index = z.getSparseGrid()->getIndexOfTable(i);

                MultiDimFiveCompass mc_outer;
                bool todo[mc_outer.getMaxShift()];

                IndexDimension index_mc[mc_outer.getMaxShift()];
                double basis_coeff_mc[mc_outer.getMaxShift()];

                for (int i = 0; i < mc_outer.getMaxShift(); i++) {
                    todo[i] = false;
                }
                Depth T(Index);
                stencilClass.initialize(T);


                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                    double basis_coeff = 1.0;
                    bool next = true;
                    IndexDimension J;
                    J = Index.nextFive(&mc, T, &basis_coeff, &next);


                    if (next) {
                        todo[mc.getShiftNumber()] = true;
                        index_mc[mc.getShiftNumber()] = J;
                        basis_coeff_mc[mc.getShiftNumber()] = basis_coeff;

                        v.setValue(J, basis_coeff);


                    }


                }

                double coeff = 0.0;

                #pragma omp parallel for schedule(dynamic)
                for (int j = 0; j < mc_outer.getMaxShift(); j++) {


                    if (todo[j]) {

                        IndexDimension J = index_mc[j];
                        double stencil = applyStencilLocal(J, v, T, stencilClass);


#pragma omp critical
                        stencil *= basis_coeff_mc[j];

#pragma omp critical
                        coeff += stencil;


                    }

                }


                v = 0.0;

                z.setValue(i, coeff);
#ifdef MY_MPI_ON
                }
#endif
            }


        }

#ifdef MY_MPI_ON
        MPI_Allreduce(MPI_IN_PLACE,z.getDatatableVector(),z.getSparseGrid()->getMaximalOccupiedSecondTable(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

    }


    void precond2(StencilClass stencilClass){

#pragma omp parallel
        {
            VectorSparseG v(*z.getSparseGrid());
            #pragma omp for schedule(runtime)
            for (size_t i = 0; i < z.getSparseGrid()->getMaximalOccupiedSecondTable(); i++) {
                if (z.getSparseGrid()->getActiveTable()[i]) {

                    IndexDimension Index = z.getSparseGrid()->getIndexOfTable(i);

                    MultiDimFiveCompass mc_outer;
                    bool todo[mc_outer.getMaxShift()];


                    IndexDimension index_mc[mc_outer.getMaxShift()];
                    double basis_coeff_mc[mc_outer.getMaxShift()];

                    for (int j = 0; j < mc_outer.getMaxShift(); j++) {
                        todo[j] = false;
                    }
                    Depth T(Index);


                    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                        double basis_coeff = 1.0;
                        bool next = true;
                        IndexDimension J;
                        J = Index.nextFive(&mc, T, &basis_coeff, &next);


                        if (next) {

                            todo[mc.getShiftNumber()] = true;
                            index_mc[mc.getShiftNumber()] = J;
                            basis_coeff_mc[mc.getShiftNumber()] = basis_coeff;

                            v.setValue(J, basis_coeff);


                        }


                    }

                    double coeff = 0.0;


                    for (int j = 0; j < mc_outer.getMaxShift(); j++) {


                        if (todo[j]) {

                            IndexDimension J = index_mc[j];
                            double stencil = applyStencilLocal(J, v, T, stencilClass);


                            stencil *= basis_coeff_mc[j];


                            coeff += stencil;



                        }

                    }


                    v = 0.0;

                    z.setValue(i, coeff);

                }
            }


        }





    }


    void precond3(StencilClass stencilClass){
#ifdef MY_MPI_ON
        int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif



        int iter_i=0;
        DepthList depthList(*z.getSparseGrid());




        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {


            Depth T = *it;

            SingleDepthHashGrid &depthGrid = z.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
            const auto &mapping = depthGrid._mapPosToGridPos;
            const auto mapping_size = mapping.size();


            if (depthGrid.getNumberOfEntries() > 0) {


#ifdef MY_MPI_ON
                if (world_rank == iter_i  % world_size ) {
#endif

                stencilClass.initialize(T);
#pragma omp parallel
                {
                    VectorSparseG v(*z.getSparseGrid());
                    #pragma omp for schedule(dynamic)
                    for (size_t i = 0; i < mapping_size; i++) {



                        if (z.getSparseGrid()->getActiveTable()[mapping[i]]) {
                            IndexDimension Index = depthGrid._map.getIndexOfTable(i);

                            MultiDimFiveCompass mc_outer;
                            bool todo[mc_outer.getMaxShift()];

                            IndexDimension index_mc[mc_outer.getMaxShift()];
                            double basis_coeff_mc[mc_outer.getMaxShift()];

                            for (int j = 0; j < mc_outer.getMaxShift(); j++) {
                                todo[j] = false;
                            }


                            for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                                double basis_coeff = 1.0;
                                bool next = true;
                                IndexDimension J;
                                if (boundary) {
                                    J = Index.nextFive_Neumann(&mc, T, &basis_coeff, &next);
                                } else
                                    J = Index.nextFive(&mc, T, &basis_coeff, &next);


                                if (next) {
                                    todo[mc.getShiftNumber()] = true;
                                    index_mc[mc.getShiftNumber()] = J;
                                    basis_coeff_mc[mc.getShiftNumber()] = basis_coeff;

                                    v.setValue(J, basis_coeff);


                                }


                            }

                            double coeff = 0.0;


                            for (int j = 0; j < mc_outer.getMaxShift(); j++) {


                                if (todo[j]) {

                                    IndexDimension J = index_mc[j];
                                    double stencil = applyStencilLocal(J, v, T, stencilClass);
                                    stencil *= basis_coeff_mc[j];


                                    coeff += stencil;


                                }

                            }


                            v = 0.0;

                            z.setValue(mapping[i], coeff);

                        }



                    }


                }


#ifdef MY_MPI_ON
                }
#endif
            }

            iter_i++;
        }

#ifdef MY_MPI_ON

        MPI_Allreduce(MPI_IN_PLACE,z.getDatatableVector(),z.getSparseGrid()->getMaximalOccupiedSecondTable(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif



    }


    void precond4(StencilClass stencilClass){
#pragma omp parallel
        {
        VectorSparseG v(z.getSparseGrid());
        VectorSparseG u(z.getSparseGrid());
#pragma omp for schedule(dynamic)
            for (unsigned long k = 0; k < z.getSparseGrid()->getMaximalOccupiedSecondTable(); k++) {

                if (z.getSparseGrid()->getActiveTable()[k]) {


                    IndexDimension Index = z.getSparseGrid()->getIndexOfTable(k);
                    Depth Tneu(Index);
                    stencilClass.initialize(Tneu);


                    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                        double basis_coeff = 1.0;
                        bool next = true;
                        IndexDimension J;
                        if (boundary) {
                            J = Index.nextFive_Neumann(&mc, Tneu, &basis_coeff, &next);
                        } else
                            J = Index.nextFive(&mc, Tneu, &basis_coeff, &next);


                        if (next) {
                            v.setValue(J, basis_coeff);
                        }




                    }


                    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                        double basis_coeff = 1.0;
                        bool next = true;
                        IndexDimension J;
                        if (boundary) {
                            J = Index.nextFive_Neumann(&mc, Tneu, &basis_coeff, &next);
                        } else
                            J = Index.nextFive(&mc, Tneu, &basis_coeff, &next);


                        if (next) {
                            basis_coeff = 0.0;

                            // basis_coeff = CalcStencilValue_Boundary(J, Tneu, v, Stiffness, stencil);
                            basis_coeff = applyStencilLocal(J,v,Tneu,stencilClass);

                            u.setValue(J, basis_coeff);
                        }

                    }


                    double coeff = 0.0;
                    for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                        double basis_coeff = 1.0;
                        bool next = true;
                        IndexDimension J;
                        if (boundary) {
                            J = Index.nextFive_Neumann(&mc, Tneu, &basis_coeff, &next);
                        } else
                            J = Index.nextFive(&mc, Tneu, &basis_coeff, &next);

                        if (next) {
                            coeff = coeff + u.getValue(J) * basis_coeff;

                        }
                        v.setValue(J, 0.0);
                        u.setValue(J, 0.0);

                    }
                    v=0.0;
                    u=0.0;

                    z.setValue(k, coeff);

                }
            }

    }

    }



    void precondLocalStiffnessMatrix(StencilClass stencilClass){
        int world_rank;
#ifdef MY_MPI_ON

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif



        int iter_i=0;
        DepthList depthList(*z.getSparseGrid());




        for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {


            Depth T = *it;
            int node = stencilClass.getNodeForDepth(T);
            if(world_rank==node){

                SingleDepthHashGrid &depthGrid = z.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
                const auto &mapping = depthGrid._mapPosToGridPos;
                const auto mapping_size = mapping.size();


                if (depthGrid.getNumberOfEntries() > 0) {


#pragma omp parallel
                    {
                        VectorSparseG v(*z.getSparseGrid());
                        VectorSparseG v_output(*z.getSparseGrid());

#pragma omp for schedule(dynamic)
                        for (size_t i = 0; i < mapping_size; i++) {


                            if (z.getSparseGrid()->getActiveTable()[mapping[i]]) {
                                IndexDimension Index = depthGrid._map.getIndexOfTable(i);


                                MultiDimFiveCompass mc_outer;
                                bool todo[mc_outer.getMaxShift()];

                                IndexDimension index_mc[mc_outer.getMaxShift()];
                                double basis_coeff_mc[mc_outer.getMaxShift()];

                                for (int j = 0; j < mc_outer.getMaxShift(); j++) {
                                    todo[j] = false;
                                }

                                std::vector<IndexDimension> IndicesVector;
                                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                                    double basis_coeff = 1.0;
                                    bool next = true;
                                    IndexDimension J;
                                    if (boundary) {
                                        J = Index.nextFive_Neumann(&mc, T, &basis_coeff, &next);
                                    } else
                                        J = Index.nextFive(&mc, T, &basis_coeff, &next);


                                    if (next) {
                                        todo[mc.getShiftNumber()] = true;
                                        index_mc[mc.getShiftNumber()] = J;
                                        basis_coeff_mc[mc.getShiftNumber()] = basis_coeff;

                                        v.setValue(J, basis_coeff);
                                        IndicesVector.push_back(J);


                                    }


                                }

                                stencilClass.applyLocalStiffnessMatricesOnIndices_onNode(v, v_output, T, IndicesVector);

                                double coeff = 0.0;


                                for (int j = 0; j < mc_outer.getMaxShift(); j++) {


                                    if (todo[j]) {

                                        IndexDimension J = index_mc[j];
                                        double stencil;


                                        stencil = v_output.getValue(J);


                                        stencil *= basis_coeff_mc[j];


                                        coeff += stencil;


                                    }

                                }


                                v = 0.0;
                                v_output = 0.0;


                                z.setValue(mapping[i], coeff);

                            }



                    }


                    }


                }

            }

            iter_i++;
        }

#ifdef MY_MPI_ON

        MPI_Allreduce(MPI_IN_PLACE,z.getDatatableVector(),z.getSparseGrid()->getMaximalOccupiedSecondTable(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

    }










    double applyStencilLocal(IndexDimension &Index, const VectorSparseG &input,Depth &T, StencilClass stencilClass) {
        double val = 0.0;
        IndexDimension NextIndex;
        //stencilClass.initialize(T);
        for (MultiDimCompass mc; mc.goon(); ++mc) {
            NextIndex = Index.nextThree2(&mc, T.returnTiefen());

            if(mc.goon()) {
                double value = 0.0;
                value = stencilClass.returnValue(Index, mc);

                val = val + value * input.getValue(NextIndex);


            }

        }
        return val;
    };



    void apply(VectorSparseG* u){
        *u = ((*u)*z);
    }

    void apply_inverse(VectorSparseG* u){
        for (unsigned long k = 0; k < z.getSparseGrid()->getMaximalOccupiedSecondTable(); k++) {
            if(z.getSparseGrid()->getActiveTable()[k]) {
                double val1 = u->getValue(k);
                double val2 = z.getValue(k);
                u->setValue(k,double(val1/val2));


            }
        }
    }


    void apply_inverse_shift(VectorSparseG* u,VectorSparseG& z_right,double lambda){
        for (unsigned long k = 0; k < z.getSparseGrid()->getMaximalOccupiedSecondTable(); k++) {
            if(z.getSparseGrid()->getActiveTable()[k]) {
                double val1 = u->getValue(k);
                double val2 = z.getValue(k)-lambda*z_right.getValue(k);
                u->setValue(k,double(val1/val2));


            }
        }
    }

    template<class Stencil>
    void apply_inverse_shift(VectorSparseG* u,Preconditioning<Stencil>& P,double lambda){
        for (unsigned long k = 0; k < z.getSparseGrid()->getMaximalOccupiedSecondTable(); k++) {
            if(z.getSparseGrid()->getActiveTable()[k]) {
                double val1 = u->getValue(k);
                double val2 = z.getValue(k)-lambda*P.getValue(k);
                u->setValue(k,double(val1/val2));


            }
        }
    }

    void apply_NEU(VectorSparseG* u){
        for (unsigned long k = 0; k < z.getSparseGrid()->getMaximalOccupiedSecondTable(); k++) {
            if(z.getSparseGrid()->getActiveTable()[k]) {
                double val1 = u->getValue(k);
                double val2 = z.getValue(k);
                u->setValue(k, double(val1/ sqrt(val2)));
            }
        }
    }

    double getValue(unsigned long k){
        return z.getValue(k);
    }

    VectorSparseG getZ(){
        return z;
    }


private:
    VectorSparseG& z;
    bool boundary;



};





class CG {
public:
    static bool solveNeumann(int maxIteration,
                             double eps,
                             VectorSparseG &x,

                             VectorSparseG &f,



                             int *iterations,
                             MatrixVectorInhomogen& matrix);

    template<class Stencil>
    static bool solveHomogenP(double eps,
                             VectorSparseG &x,
                             VectorSparseG &f,
                             int *iterations,
                             VectorSparseG &usol,
                             MatrixVectorHomogen& matrix,Stencil stencil,double* time_precon,Preconditioning<Stencil>& P);

    template<class StencilLeft, class StencilRight>
 static bool solveHomogenShift(double eps,
                          VectorSparseG &x,
                          VectorSparseG &f,
                          int *iterations,
                          VectorSparseG &usol,
                          MatrixVectorHomogen& matrix,StencilLeft stencil_left,StencilRight stencil_right,double* time_precon, double lambda,Preconditioning<StencilLeft>& Pleft,  Preconditioning<StencilRight>& Pright);


    template<class StencilLeft, class StencilRight>
 static bool solveHomogenShift2(double eps,
                          VectorSparseG &x,
                          VectorSparseG &f,
                          int *iterations,
                          VectorSparseG &usol,
                          MatrixVectorHomogen& matrix,StencilLeft stencil_left,StencilRight stencil_right,double* time_precon, double lambda,Preconditioning<StencilLeft>& Pleft,  Preconditioning<StencilRight>& Pright);


    template<class Stencil>
static bool solveHomogen(double eps,
                         VectorSparseG &x,
                         VectorSparseG &f,
                         int *iterations,
                         VectorSparseG &usol,
                         MatrixVectorHomogen& matrix,Stencil stencil,double* time_precon);


    static const int maxIteration = 150;
    static double error[MaximumDepth][maxIteration];
    static int count;


    static void printError() {
        ofstream Datei;
        Datei.open("../results/error_CG.gnu", std::ios::out);

        for (int i = 0; i < MaximumDepth; i++) {
            if (error[i][0] != 0) {
                Datei << "\t" << error[i][0];
                cout << "\t" << error[i][0];
            }
        }
        Datei << endl;
        cout << endl;
        for (int i = 1; i < maxIteration; i++) {
            Datei << i;
            cout << i;
            for (int j = 0; j < MaximumDepth; j++) {
                if (error[j][0] != 0) {
                    Datei << "\t" << error[j][i];
                    cout << "\t " << error[j][i];
                }


            }

            Datei << "\n";
            cout << "\n";

        }
        Datei << endl;
        Datei.close();
    }
};


class OperatorT {
public:
    static void apply(VectorSparseG &input, VectorSparseG &output);

    static void applyTranspose(VectorSparseG &input, VectorSparseG &output);
};


//INLINE TEMPLATE FUNCTIONS
template<class Stencil>
bool
CG::solveHomogenP(double eps, VectorSparseG &x,VectorSparseG &f,int *iterations,
                 VectorSparseG &usol, MatrixVectorHomogen& matrix,Stencil stencil,double* time_precon, Preconditioning<Stencil>& P) {

        // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);

    double tau;
    double delta;
    double delta_prime=1.0;
    double beta;


    AdaptiveSparseGrid* grid = (x.getSparseGrid());
    //ListOfDepthOrderedSubgrids list(*grid);

    if(!(x.getSparseGrid()->getKey()==f.getSparseGrid()->getKey()))exit(1);



    VectorSparseG z(grid);
    VectorSparseG r(grid);
    VectorSparseG g(grid);
    VectorSparseG d(grid);






    VectorSparseG h(grid);
    double j = 1.0;

    // solves  Ax = f;


    // r = A*x - f;

    // Start measuring time
    struct timeval begin_all, end_all;
    gettimeofday(&begin_all, 0);





    //dirichlet_grid.completeDirichletGrid();
    //Preconditioning<Stencil> P(z,stencil);




    matrix.multiplication<Stencil>(x, r, stencil);



    // r = r - f = Ax - f;
    r = f - r;

    h = r;



    P.apply_inverse(&h);



    d = -1.0*h;




    delta = product(r, h);




    int k = 0;

    for (int i = 1; i <= maxIteration && delta > eps; ++i) {
        grid->WorkOnHangingNodes = false;
        k = i;

        g = 0.0;


        matrix.multiplication<Stencil>(d, g,stencil);


        tau = double(delta / product(d, g));

        r = r + (tau * g);
        x = x - (tau * d);

        h = j * r;

        P.apply_inverse(&h);

        delta_prime = product(r, h);
        if (delta_prime < eps) break;


        beta = delta_prime / delta;
        delta = delta_prime;
        d = -1.0*h + beta * d;


    }


    //cout << "CG finished after " << k << " iterations and residuum delta = " << sqrt(delta_prime) << endl;
    *iterations = k;
    count++;

    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double duration_def = seconds + microseconds*1e-6;
    *time_precon = duration_def;
    return true;
}



//INLINE TEMPLATE FUNCTIONS
template<class StencilLeft,class StencilRight>
bool
CG::solveHomogenShift(double eps, VectorSparseG &x,VectorSparseG &f,int *iterations,
                 VectorSparseG &usol, MatrixVectorHomogen& matrix,StencilLeft stencil_left,StencilRight stencil_right, double* time_precon,double lambda,Preconditioning<StencilLeft>& Pleft,  Preconditioning<StencilRight>& Pright) {

    // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);


    double tau;
    double delta;
    double delta_prime=1.0;
    double beta;


    AdaptiveSparseGrid* grid = (x.getSparseGrid());


    if(!(x.getSparseGrid()->getKey()==f.getSparseGrid()->getKey()))exit(1);



    VectorSparseG z_right(grid);
    z_right = Pright.getZ();
    VectorSparseG r(grid);
    VectorSparseG r_right(grid);
    VectorSparseG g(grid);
    VectorSparseG g_right(grid);
    VectorSparseG d(grid);






    VectorSparseG h(grid);
    double j = 1.0;

    // solves  (A-lambda*M)x = f;


    // r = A*x - f;











    matrix.multiplication<StencilLeft>(x, r, stencil_left);
    matrix.multiplication<StencilRight>(x,r_right,stencil_right);
    r = r- lambda*r_right;


    // r = r - f = Ax - f;
    r = f - r;

    h = r;



   // Pleft.apply_inverse_shift(&h,z_right,lambda);



    d = -1.0*h;




    delta = product(r, h);




    int k = 0;

    for (int i = 1; i <= maxIteration && delta > eps; ++i) {
        grid->WorkOnHangingNodes = false;
        k = i;

        g = 0.0;



        matrix.multiplication<StencilLeft>(d, g,stencil_left);
        matrix.multiplication<StencilRight>(d, g_right,stencil_right);
        g = g - lambda*g_right;


        tau = double(delta / product(d, g));

        r = r + (tau * g);
        x = x - (tau * d);

        h = j * r;

        //Pleft.apply_inverse_shift(&h,z_right,lambda);

        delta_prime = product(r, h);
        if (delta_prime < eps) break;


        beta = delta_prime / delta;
        delta = delta_prime;
        d = -1.0*h + beta * d;


    }


    //cout << "CG finished after " << k << " iterations and residuum delta = " << sqrt(delta_prime) << endl;
    *iterations = k;
    count++;

    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double duration_def = seconds + microseconds*1e-6;
    *time_precon = duration_def;
    return true;
}


template<class StencilLeft,class StencilRight>
bool
CG::solveHomogenShift2(double eps, VectorSparseG &x,VectorSparseG &f,int *iterations,
                 VectorSparseG &usol, MatrixVectorHomogen& matrix,StencilLeft stencil_left,StencilRight stencil_right, double* time_precon,double lambda,Preconditioning<StencilLeft>& Pleft,  Preconditioning<StencilRight>& Pright) {


    AdaptiveSparseGrid* grid = (x.getSparseGrid());
    VectorSparseG r(grid),Ax(grid),Bx(grid),z(grid),p(grid),Sp(grid);

    matrix.multiplication<StencilLeft>(x, Ax, stencil_left);
    matrix.multiplication<StencilRight>(x, Bx, stencil_right);

    r = f - (Ax-lambda*Bx);
    z = r;
    Pleft.apply_inverse_shift(&z,Pright,lambda);
    p = z;

    double rs_old = product(r,p);

    for (int i = 1; i <= maxIteration; ++i)
    {
        matrix.multiplication<StencilLeft>(p, Ax, stencil_left);
        matrix.multiplication<StencilRight>(p, Bx, stencil_right);
        Sp = Ax - lambda*Bx;
        double alpha = rs_old/product(p,Sp);

        x = x + alpha*p;
        r = r - alpha*Sp;

        double res_norm = sqrt(product(r,r));
        if (res_norm < eps) break;

        z=r;
        Pleft.apply_inverse_shift(&z,Pright,lambda);

        double rs_new = product(r,z);
        double beta = rs_new/rs_old;

        p = z + beta*p;
        rs_old = rs_new;
        *iterations = i;

    }
    return true;
}



//INLINE TEMPLATE FUNCTIONS
template<class Stencil>
bool
CG::solveHomogen(double eps, VectorSparseG &x,VectorSparseG &f,int *iterations,
                 VectorSparseG &usol, MatrixVectorHomogen& matrix,Stencil stencil,double* time_precon) {

        // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);

    double tau;
    double delta;
    double delta_prime=1.0;
    double beta;


    AdaptiveSparseGrid* grid = (x.getSparseGrid());
    ListOfDepthOrderedSubgrids list(*grid);

    if(!(x.getSparseGrid()->getKey()==f.getSparseGrid()->getKey()))exit(1);



    VectorSparseG z(grid);
    VectorSparseG r(grid);
    VectorSparseG g(grid);
    VectorSparseG d(grid);




    VectorSparseG h(grid);
    double j = 1.0;

    // solves  Ax = f;


    // r = A*x - f;

    // Start measuring time
    struct timeval begin_all, end_all;
    gettimeofday(&begin_all, 0);





    //dirichlet_grid.completeDirichletGrid();
    Preconditioning<Stencil> P(z,stencil);




    matrix.multiplication<Stencil>(x, r, stencil);



    // r = r - f = Ax - f;
    r = f - r;

    h = r;



    P.apply_inverse(&h);


    d = -1.0*h;




    delta = product(r, h);




    int k = 0;

    for (int i = 1; i <= maxIteration && delta > eps; ++i) {
        grid->WorkOnHangingNodes = false;
        k = i;

        g = 0.0;


        matrix.multiplication<Stencil>(d, g,stencil);


        tau = double(delta / product(d, g));

        r = r + (tau * g);
        x = x - (tau * d);

        h = j * r;

        P.apply_inverse(&h);
        delta_prime = product(r, h);
        if (delta_prime < eps) break;


        beta = delta_prime / delta;
        delta = delta_prime;
        d = -1.0*h + beta * d;


    }


    //cout << "CG finished after " << k << " iterations and residuum delta = " << sqrt(delta_prime) << endl;
    *iterations = k;
    count++;
    
    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double duration_def = seconds + microseconds*1e-6;
    *time_precon = duration_def;
    return true;
}

#endif //GRUN_CG_METHOD_H
