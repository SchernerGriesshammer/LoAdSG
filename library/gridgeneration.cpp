#include <iomanip>

#include "source/sparseEXPDE.h"


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

double var_coeff( double* x)
{
    const double sum = (x[0]*x[0])+(x[1]*x[1])+(x[2]*x[2]);
    return 1.0/sum;
}



int main(int argc, char **argv) {

#ifdef MY_MPI_ON
    int rank, num_tasks;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
#endif


    AdaptiveSparseGrid regular_grid;
    AdaptiveSparseGrid grid2;

    cout << DimensionSparseGrid    <<endl;

    int d = 6;


    IndexDimension centerPoint;

    regular_grid.AddPoint(centerPoint.leftSon(0).leftSon(0));
    regular_grid.completeDirichletGrid();
    regular_grid.PrintActiveHanging(4);

    VectorSparseG u(regular_grid),v(regular_grid);
    u=1.0;
    calcNodalByPrew(u,v);
    u=0.0;
    v.PrintDouble(3);
    calcPrewByNodal(u,v);

    u.PrintDouble(3);


    return 0;
    grid2.AddRecursiveSonsOfPoint(centerPoint,6);

    int size=grid2.getMaximalOccupiedSecondTable();

    MultiLevelAdaptiveSparseGrid mgrid(&regular_grid);
    MatrixVectorHomogen m(regular_grid,mgrid,true,false);
    Poisson stencilPoisson(regular_grid);


    VariableCoefficientMass helm_holtz(regular_grid);

    double matrix[size][size];
    VectorSparseG x(regular_grid);
    VectorSparseG b(regular_grid);

    calcPrewByNodal(x,b);

    for (unsigned long k=0; k<grid2.getMaximalOccupiedSecondTable();k++)
    {
        IndexDimension p=grid2.getIndexOfTable(k);
        x.setValue(p,1.0);

        m.multiplication(x,b,helm_holtz);
        for (unsigned long j=0;j< grid2.getMaximalOccupiedSecondTable();j++)
        {
            IndexDimension q=grid2.getIndexOfTable(j);

            matrix[k][j]=b.getValue(q);
        }

        x=0.0;
        b=0.0;

    }



    std::ofstream datei("matrix_data_helmsparse7.txt");

    if (datei.is_open()) {
        // Hohe Präzision sicherstellen, damit double-Werte nicht gerundet werden
        datei << std::setprecision(15);

        for (int i = 0; i < grid2.getMaximalOccupiedSecondTable(); ++i) {
            for (int j = 0; j < grid2.getMaximalOccupiedSecondTable(); ++j) {
                datei << matrix[i][j];
                if (j < grid2.getMaximalOccupiedSecondTable() - 1) datei << " "; // Leerzeichen zwischen Spalten
            }
            datei << "\n"; // Zeilenumbruch nach jeder Zeile
        }
        datei.close();
        std::cout << "Matrix erfolgreich gespeichert." << std::endl;
    }

    






#ifdef MY_MPI_ON
    MPI_Finalize();
#endif
    return 0;
}