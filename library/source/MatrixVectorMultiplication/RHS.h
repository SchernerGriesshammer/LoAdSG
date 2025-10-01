//
// Created by to35jepo on 3/13/23.
//

#ifndef SGRUN_RHS_H
#define SGRUN_RHS_H
#include "../extemp/vector.h"
#include "../extemp/multilevelvector.h"
#include "../stencils/MassStencil.h"
#include "MatrixVectorInhomogen.h"
#include "../iterator/RectangularIterator.h"
#include "../stencils/InterfaceIntegration/IntegratorRHS/integratorRHS.h"



void CalcHierarchicalBasisForRHS(VectorSparseG &u);


class InHomoBoundaryRHS{
public:
    InHomoBoundaryRHS(AdaptiveSparseGrid& grid, MultiLevelAdaptiveSparseGrid& mgrid): depthList(grid), nodalMultiLevel(mgrid), nodal(grid){
        bool prolongations[DimensionSparseGrid]={0};
        depthList.sortDepthsNeumann(prolongations);
    };

    void dualTransform(VectorSparseG& neumann, VectorSparseG& dirichlet){
        dirichlet=0.0;
        IndexDimension minI;
        IndexDimension maxI;

        std::list<Depth>* sortedDepths = depthList.getSortierteTiefen();
        for(Depth T: *sortedDepths) nodalMultiLevel.setMultiLevelValues2(neumann, T);


        for(int d=0; d < DimensionSparseGrid; d++) {

            for (Depth T: *sortedDepths) {


                Depth T_fine=T;
                int t = T.at(d) + 1;
                T_fine.set(t,d);


                if(depthList.isIncluded(T_fine)) {

                    Matrix M(t);

                    IndexDimension I;

                    for (int dd = 0; dd < DimensionSparseGrid; dd++) {
                        minI.replace(dd, 0);
                        maxI.replace(dd, 1);
                    }
                    maxI.replace(d,0);

                    double PNeumann[M.size_all];
                    double PNodal[M.size_all];
                    unsigned long indices[M.size_all];


                    for (RectangularIterator iter(minI, maxI, T); iter.goon();++iter) {
                        I = iter.getIndex();

                        for (int i = 0; i < M.size_all; i++) {


                            unsigned long k;
                            nodalMultiLevel.getSparseGrid()->occupied(k, I, T_fine);
                            indices[i]=k;
                            if(i%2==0){
                                nodalMultiLevel.getSparseGrid()->occupied(k, I, T);
                                PNeumann[i]=nodalMultiLevel.getValue(k);
                            } else {
                                PNeumann[i]=nodalMultiLevel.getValue(indices[i]);
                            }

                            I = I.nextRight(d, T_fine.at(d));
                        }

                        M.solve(PNodal,PNeumann);

                        for (int i = 0; i < M.size_all; i++) {

                            nodalMultiLevel.setValue(indices[i], PNodal[i]);

                        }

                    }

                }

            }

        }

        for(Depth T: *sortedDepths){
            if(T>>0) {
                nodal.setMultiLevelValues(nodalMultiLevel, T);
                ConvertDualNodalToPrewavelet(nodal, dirichlet, T);

            }
        }
    }

    class Matrix{
    public:
        Matrix(int t){


            size_all = POW2(t)+1;
            if (size_all < 6 && size_all >= 1) {
                size = size_all*size_all;
            }else {
                size = size_all * 5;
            }

            M = new double[size];
            L = new double[size];
            R = new double[size];


            for(int k=0; k<size; k++){
                M[k]=0.0;
                L[k]=0.0;
                R[k]=0.0;
            }


            if(t>1) {
                for (int k = 2; k < size_all - 1; k += 2) {
                    M[convert(k, k, size_all)] = 1.0;
                    M[convert(k, k - 1, size_all)] = 0.5;
                    M[convert(k, k + 1, size_all)] = 0.5;
                }
                for (int k = 3; k < size_all - 2; k += 2) {
                    M[convert(k, k, size_all)] = 1.0;
                    M[convert(k, k - 1, size_all)] = -0.6;
                    M[convert(k, k - 2, size_all)] = 0.1;
                    M[convert(k, k + 1, size_all)] = -0.6;
                    M[convert(k, k + 2, size_all)] = 0.1;
                }

                M[convert(0, 0, size_all)] = 1.0;
                M[convert(0, 1, size_all)] = 0.5;
                M[convert(1, 0, size_all)] = -1.2;
                M[convert(1, 1, size_all)] = 1.1;
                M[convert(1, 2, size_all)] = -0.6;
                M[convert(1, 3, size_all)] = 0.1;

                M[convert(size_all - 1, size_all - 1, size_all)] = 1.0;
                M[convert(size_all - 1, size_all - 2, size_all)] = 0.5;
                M[convert(size_all - 2, size_all - 1, size_all)] = -1.2;
                M[convert(size_all - 2, size_all - 2, size_all)] = 1.1;
                M[convert(size_all - 2, size_all - 3, size_all)] = -0.6;
                M[convert(size_all - 2, size_all - 4, size_all)] = 0.1;
            }
            if(t==1){
                M[convert(0,0,size_all)]=1.0;
                M[convert(0,1,size_all)]=0.5;
                M[convert(1,0,size_all)]=-1.0;
                M[convert(1,1,size_all)]=1.0;
                M[convert(1,2,size_all)]=-1.0;
                M[convert(2,1,size_all)]=0.5;
                M[convert(2,2,size_all)]=1.0;
            }

            if(t==0){  M[convert(0,0,size_all)]=1.0;  M[convert(1,1,size_all)]=1.0;}

            LRcreated = false;

        }
        ~Matrix(){

            delete[] M;
            delete[] L;
            delete[] R;

        }
        void print() {
            int I = size_all;
            for (int i = 0; i < I; i++) {
                for (int j = 0; j < I; j++) {
                    if (i <= j + 2 && j <= i + 2) {
                        cout << "\t" << M[convert(i, j, I)];
                    } else { cout << "\t"; }
                }
                cout << "\n";
            }
            cout << endl;

        };

        void solve(double *x, double *b) {
            if (size_all < 1) {
                x[0] = b[0];
                return;
            }
            if (!LRcreated) createLR();
            double y[size_all];
            for (int i = 0; i < size_all; i++) {
                y[i] = b[i];
            }

            forward(y, b, L, size_all);
            backward(x, y, R, size_all);
        };


        int size_all;
    private:
        void createLR() {
            for (int i = 0; i < size_all; i++) {
                for (int j = 0; j < size_all; j++) {
                    if (abs(i - j) < 3) {
                        R[convert(i, j, size_all)] = M[convert(i, j, size_all)];
                    }
                }
            }

            LR(L, R, size_all);
            LRcreated = true;
        };

        int convert(int i,int j,int I){


            int savei = i;
            if (I < 6) { return i * I + j; }
            if (j > 2 && j < I - 2) {
                i = i - ((j - 3) + 1);
            }
            if (j == I - 2) i = i - j + 3;
            if (j == I - 1) i = i - j + 4;

            if (I > 5) {
                if (i > 4 || i < 0) {
                    cout << "error" << endl;
                    cout << " i " << i << " j  " << j << endl;
                    cout << "save i  " << savei << "   I " << I << endl;
                }
            }


            return i * I + j;
        }


        double *M;
        double *R;
        double *L;
        bool LRcreated;

        int size;

    };




private:

    DepthList depthList;
    MultiLevelVector nodalMultiLevel;
    VectorSparseG nodal;



};








class InHomoBoundaryRHSHierarchical{
public:
    InHomoBoundaryRHSHierarchical(AdaptiveSparseGrid& grid):nodal(grid){};

    void multiply(VectorSparseG& hierarchical, VectorSparseG& dirichlet);

    void multiply_mass(VectorSparseG& hierarchical, VectorSparseG& dirichlet);

    void multiply_mass_coeff(VectorSparseG &hierarchical, VectorSparseG &dirichlet);

    template<class F>
    void multiply_mass_coeff_interface(VectorSparseG &hierarchical, VectorSparseG &dirichlet,const F& var_coeff) {


        DepthList depthList(*dirichlet.getSparseGrid());
        DepthList depthList1(*hierarchical.getSparseGrid());

            for (auto it = depthList.begin_all(); it != depthList.end_all(); ++it) {
                Depth T = *it;
                auto iter = GetNextSmallerDepthsIterator(T);
                do {
                    Depth Tlocal = *iter;
                    if (Tlocal >> 0) {

                        SingleDepthHashGrid &depthGrid = dirichlet.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                                Tlocal);
                        const auto &mapping = depthGrid._mapPosToGridPos;

                        if (depthGrid.getNumberOfEntries() > 0) {


                            for (size_t i = 0; i < mapping.size(); i++) {

                                IndexDimension INodal = depthGrid._map.getIndexOfTable(i);

                                //double sum = 0.0;
                                double sum_omp = 0.0;
                                double sum = 0.0;


                                for (auto it_inner = depthList1.begin_all();
                                     it_inner != depthList1.end_all(); ++it_inner) {

                                    Depth THier = *it_inner;
                                    if (!(THier >> 0)) {

                                        SingleDepthHashGrid &depthGrid_hier = hierarchical.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(
                                                THier);
                                        const auto &mapping_hier = depthGrid_hier._mapPosToGridPos;
                                        const auto &map_hier = depthGrid_hier._map;
                                        auto end = mapping_hier.size();


                                        if (depthGrid_hier.getNumberOfEntries() > 0){
                                            #pragma omp parallel reduction(+:sum_omp)
                                            {

                                                #pragma omp for schedule(runtime)
                                                for (size_t j = 0; j < end; j++) {
                                                    //IndexDimension IHier = depthGrid_hier._map.getIndexOfTable(j);
                                                    IndexDimension IHier = map_hier.getIndexOfTable(j);

                                                    array<double, DimensionSparseGrid> p_u{}, p_v{}, l_b{}, r_b{}, h_u{}, h_v{};
                                                    array<int, DimensionSparseGrid> d_u{}, d_v{};

                                                    for (int d = 0; d < DimensionSparseGrid; d++) {
                                                        p_u[d] = INodal.coordinate(d);
                                                        p_v[d] = IHier.coordinate(d);
                                                        h_u[d] = 1.0 / double(POW2(T.at(d)));
                                                        h_v[d] = 1.0 / double(POW2(THier.at(d)));
                                                        l_b[d] = max(p_u[d] - h_u[d], p_v[d] - h_v[d]);
                                                        l_b[d] = max(l_b[d], 0.0);
                                                        r_b[d] = min(p_u[d] + h_u[d], p_v[d] + h_v[d]);
                                                        r_b[d] = min(r_b[d], 1.0);
                                                        d_u[d] = T.at(d);
                                                        d_v[d] = THier.at(d);
                                                    }

                                                    IntegratorRHS<const F, DimensionSparseGrid> integratorRhs(var_coeff,
                                                             1e-4, 0, 3, 1.0e8);
                                                    double val = integratorRhs.integration(l_b, r_b,
                                                                                           p_u, p_v,
                                                                                           h_u,
                                                                                           h_v,
                                                                                           d_u, d_v);

                                                    //#pragma omp critical
                                                    sum_omp += val *
                                                               hierarchical.getValue(mapping_hier[j]);


                                                }

                                            }

                                        }
                                    }
                                }

                                              /*  if(abs(sum-sum_omp)>1e-10){
                                                    cout << sum << " " << sum_omp << "  " << sum_omp-sum << endl;
                                                    exit(1);

                                                }*/
                                nodal.setValue(mapping[i], sum_omp);

                            }
                        }
                    }


                } while (iter.next());

                ConvertToPrewavelet2(nodal, dirichlet, T);
                nodal = 0.0;
            }

    }
private:
    inline void ConvertToPrewavelet2(VectorSparseG &Ax_hier, VectorSparseG &Ax, Depth &T) {
        SingleDepthHashGrid& depthGrid = Ax_hier.getSparseGrid()->getMultiDepthHashGrid()->getGridForDepth(T);
        const auto& mapping = depthGrid._mapPosToGridPos;
        // cout << mapping.size() << ", " << depthGrid.getNumberOfEntries() <<endl;
        for (size_t i = 0; i < mapping.size(); i++)
        {
            if(Ax_hier.getSparseGrid()->getActiveTable()[mapping[i]]) {
                IndexDimension I = depthGrid._map.getIndexOfTable(i);
                double coeff = 0.0;
                IndexDimension J;

                for (MultiDimFiveCompass mc; mc.goon(); ++mc) {
                    double basis_coeff = 1.0;
                    J = I.nextFiveP(&mc, T, &basis_coeff);
                    if (mc.goon())
                        coeff = coeff + Ax_hier.getValue(J) * basis_coeff;
                }

                Ax.setValue(mapping[i], coeff);
            }
        }
    }


    VectorSparseG nodal;



};

#endif //SGRUN_RHS_H
