#include <filesystem>
#include <iomanip>

#include "../source/sparseEXPDE.h"


class NodalTransformation
{
public:
    NodalTransformation(VectorSparseG& u_):u(u_){};
    void apply(VectorSparseG &output,Depth T, int d)
    {
        for (unsigned long i=0; i<u.getSparseGrid()->getMaximalOccupiedSecondTable(); i++)
        {
            IndexDimension I = u.getSparseGrid()->getIndexOfTable(i);
            bool doit=false;


            Depth TI(I);


            for (int s=0; s<DimensionSparseGrid;s++)
            {
                if (nodal[s]&&TI.at(s)<=T.at(s)) doit = true;

                else if (!nodal[s]&& TI.at(s)==T.at(s)) doit = true;

                else
                {
                    doit = false;
                    break;
                }


            }




            if (doit){
                I.PrintCoord();
                if (TI.at(d)==1) output.addToValue(i,u.getValue(i));
                else if (I.isLinksRandNah(d))
                {

                    IndexDimension center=I;
                    IndexDimension right=I.nextRight(d,T.at(d));
                    IndexDimension rightright=I.nextRight(d,T.at(d)).nextRight(d,T.at(d));

                    double val = u.getValue(center);
                    output.addToValue(center,0.9*val);
                    output.addToValue(right,-0.6*val);
                    output.addToValue(rightright,0.1*val);
                }
                else if (I.isRechtsRandNah(d))
                {

                    IndexDimension center=I;
                    IndexDimension left=I.nextLeft(d,T.at(d));
                    IndexDimension leftleft=I.nextLeft(d,T.at(d)).nextLeft(d,T.at(d));


                    double val= u.getValue(i);
                    output.addToValue(center,0.9*val);

                    output.addToValue(left,-0.6*val);
                    output.addToValue(leftleft,0.1*val);
                }else
                {

                    IndexDimension center=I;
                    IndexDimension right=I.nextRight(d,T.at(d));
                    IndexDimension rightright=I.nextRight(d,T.at(d)).nextRight(d,T.at(d));

                    double val= u.getValue(i);

                    output.addToValue(center,val);

                    output.addToValue(right,-0.6*val);
                    output.addToValue(rightright,0.1*val);
                    IndexDimension left=I.nextLeft(d,T.at(d));
                    IndexDimension leftleft=I.nextLeft(d,T.at(d)).nextLeft(d,T.at(d));


                    output.addToValue(left,-0.6*val);
                    output.addToValue(leftleft,0.1*val);

                }
            }


        }
        nodal[d]=true;
        u=output;
    }
private:
    bool nodal[DimensionSparseGrid]{};
    VectorSparseG u;

};
void CalcNodalinDir(VectorSparseG& input, VectorSparseG& output, Depth T, int d, bool* nodal)
{
    for (unsigned long i=0; i<input.getSparseGrid()->getMaximalOccupiedSecondTable(); i++)
    {
        IndexDimension I = input.getSparseGrid()->getIndexOfTable(i);
        bool doit=false;


        Depth TI(I);


        for (int s=0; s<DimensionSparseGrid;s++)
        {
            if (nodal[s]&&TI.at(s)<=T.at(s)) doit = true;

            else if (!nodal[s]&& TI.at(s)==T.at(s)) doit = true;

            else
            {
                doit = false;
                break;
            }


        }




        if (doit){
            I.PrintCoord();
            if (TI.at(d)==1) output.addToValue(i,input.getValue(i));
            else if (I.isLinksRandNah(d))
            {

                IndexDimension center=I;
                IndexDimension right=I.nextRight(d,T.at(d));
                IndexDimension rightright=I.nextRight(d,T.at(d)).nextRight(d,T.at(d));

                double val = input.getValue(center);
                output.addToValue(center,0.9*val);
                output.addToValue(right,-0.6*val);
                output.addToValue(rightright,0.1*val);
            }
            else if (I.isRechtsRandNah(d))
            {

                IndexDimension center=I;
                IndexDimension left=I.nextLeft(d,T.at(d));
                IndexDimension leftleft=I.nextLeft(d,T.at(d)).nextLeft(d,T.at(d));


                double val= input.getValue(i);
                output.addToValue(center,0.9*val);

                output.addToValue(left,-0.6*val);
                output.addToValue(leftleft,0.1*val);
            }else
            {

                IndexDimension center=I;
                IndexDimension right=I.nextRight(d,T.at(d));
                IndexDimension rightright=I.nextRight(d,T.at(d)).nextRight(d,T.at(d));

                double val= input.getValue(i);

                output.addToValue(center,val);

                output.addToValue(right,-0.6*val);
                output.addToValue(rightright,0.1*val);
                IndexDimension left=I.nextLeft(d,T.at(d));
                IndexDimension leftleft=I.nextLeft(d,T.at(d)).nextLeft(d,T.at(d));


                output.addToValue(left,-0.6*val);
                output.addToValue(leftleft,0.1*val);

            }
        }


    }
}

void tensor(Depth TX, Depth TY, AdaptiveSparseGrid& grid)
{
    cout << "tensor" << endl;
    TX.Print();
    TY.Print();
    unsigned long end = grid.getMaximalOccupiedSecondTable();
    for (unsigned long k=0; k <end; k++)
    {
        if (grid.getActiveTable()[k])
        {
            IndexDimension IX = grid.getIndexOfTable(k);

            Depth TX_lokal(IX);


            if (TX_lokal.at(0)<=TX.at(0) && TX_lokal.at(1)==TX.at(1))
            {

                for (unsigned long j=0; j < end; j++)
                {
                    IndexDimension IY = grid.getIndexOfTable(j);

                    Depth TY_lokal(IY);
                    if (TY_lokal.at(1)<=TY.at(1) &&  TY_lokal.at(0)==TY.at(0))
                    {

                        IY.replace(0,IX.getIndex(0));
                        grid.AddPoint(IY);
                    }

                }

            }
        }

    }
}


//hydrogen values
double border = 7.5;
//3e-9;
//8e-10;

double eps = 1e-10;
double cg_eps = 1e-10;


const double CUTCOEFF = 1e6;

const int mc_samples = 200;



const double varc_factor = -8.0*border;

const double rhs_factor = 8.0*border*border;

const double eig_exp = -0.5;//*rhs_factor;
const double shiftfactor =0.9;
const double shiftvalue = (shiftfactor+1.0) * (-1.0) * eig_exp*rhs_factor;


volatile double r_min = 1;
volatile double r_max = 1;

double alpha = 1.0;



double var_coeff( double* coordinates)
{
    double r = 0;
    for (int i = 0; i < DimensionSparseGrid; i++) {
        double ci = 2. * coordinates[i] - 1.; // < map from [0,1] to [-s,s]
        r += ci * ci;
    }

    r = sqrt(r);


    double result = CUTCOEFF;
    if (r!=0){
        result = min(1. / r, CUTCOEFF); // < clipping large values of 1/r around the center
    }

    result = varc_factor * result;
    //result -= 1.9*eig_exp; // + 1.9*1*4*s^2 (eig_exp==1)
    result += shiftvalue;


    // log range of sampled values
    double over_r = (r!=0) ? 1./r : 1;
    if (alpha==1 && over_r < r_min) r_min = over_r;
    if (alpha==1 && over_r > r_max) r_max = over_r;



    return alpha*result;
}



int main(int argc, char **argv) {

    int rank = 0;
    int num_tasks = 1;


#ifdef MY_MPI_ON
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#endif



    int numberMVprocesses=4;
    int numberLSprocesses=num_tasks - numberMVprocesses;


    AdaptiveSparseGrid regular_grid;


    IndexDimension centerPoint;

    regular_grid.AddRecursiveSonsOfPoint(centerPoint,5);

    VectorSparseG x(regular_grid),Ax(regular_grid),Bx(regular_grid);




    int size = regular_grid.getDOFS();

    MultiLevelAdaptiveSparseGrid mgrid(&regular_grid);
    MatrixVectorHomogen m(regular_grid,mgrid,true,true);


    double mAx[size][size];
    double diagA[size][size];
    double mBx[size][size];
    double diagB[size][size];


    Poisson poisson(regular_grid);
    HelmHoltz rhs(regular_grid);

    StencilMC<double (*)(double *)> stencilVarCoeff(regular_grid,&var_coeff,mc_samples);
    LocalStiffnessMatricesDynamicDistribution lhs(regular_grid, stencilVarCoeff,numberLSprocesses);
    lhs.addStencil(poisson);

    for (unsigned long k=0; k<regular_grid.getMaximalOccupiedSecondTable();k++)
    {
        IndexDimension p=regular_grid.getIndexOfTable(k);
        x.setValue(p,1.0);


        m.multiplication(x,Ax,lhs);
        for (unsigned long j=0;j< regular_grid.getMaximalOccupiedSecondTable();j++)
        {
            IndexDimension q=regular_grid.getIndexOfTable(j);

            mAx[k][j]=Ax.getValue(q);
        }

        diagA[k][k]=1.0/mAx[k][k];



        x=0.0;
        Ax=0.0;

    }


    for (unsigned long k=0; k<regular_grid.getMaximalOccupiedSecondTable();k++)
    {
        IndexDimension p=regular_grid.getIndexOfTable(k);
        x.setValue(p,1.0);

        m.multiplication(x,Bx,rhs);
        for (unsigned long j=0;j< regular_grid.getMaximalOccupiedSecondTable();j++)
        {
            IndexDimension q=regular_grid.getIndexOfTable(j);

            mBx[k][j]=Bx.getValue(q);
        }

        diagB[k][k]=1.0/mBx[k][k];



        x=0.0;
        Bx=0.0;

    }

if (rank==0)
{
    std::ofstream datei("mAneu5mc200cc6.txt");

    if (datei.is_open()) {
        // Hohe Präzision sicherstellen, damit double-Werte nicht gerundet werden
        datei << std::setprecision(15);

        for (int i = 0; i < regular_grid.getMaximalOccupiedSecondTable(); ++i) {
            for (int j = 0; j < regular_grid.getMaximalOccupiedSecondTable(); ++j) {
                datei << mAx[i][j];
                if (j < regular_grid.getMaximalOccupiedSecondTable() - 1) datei << " "; // Leerzeichen zwischen Spalten
            }
            datei << "\n"; // Zeilenumbruch nach jeder Zeile
        }
        datei.close();
        std::cout << "Matrix  A erfolgreich gespeichert." << std::endl;
    }



    std::ofstream datei_B("mBneu5mc200cc6.txt");

    if (datei_B.is_open()) {
        // Hohe Präzision sicherstellen, damit double-Werte nicht gerundet werden
        datei_B << std::setprecision(15);

        for (int i = 0; i < regular_grid.getMaximalOccupiedSecondTable(); ++i) {
            for (int j = 0; j < regular_grid.getMaximalOccupiedSecondTable(); ++j) {
                datei_B << mBx[i][j];
                if (j < regular_grid.getMaximalOccupiedSecondTable() - 1) datei_B << " "; // Leerzeichen zwischen Spalten
            }
            datei_B << "\n"; // Zeilenumbruch nach jeder Zeile
        }
        datei_B.close();
        std::cout << "Matrix B erfolgreich gespeichert." << std::endl;
    }
}
    






#ifdef MY_MPI_ON
    MPI_Finalize();
#endif
    return 0;
}