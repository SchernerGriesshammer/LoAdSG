#include "index.h"
#include "../sgrid/depth.h"

IndexDimension IndexDimension::nextThree_boundary(MultiDimCompass *mc, Depth &T) const {

    IndexDimension J = (*this);

    while (mc->goon()) {

        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = T.at(d);
            Richtung r = mc->getRichtung(d);
            //if(J.getDepth(d)==1) r = Mitte;



            if (r == Links) {


                if (J.isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t);
            }


            if (r == Rechts) {


                if (J.isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t);
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.
    return J;
}


/**
 * Berechnet 5^d zur Compile-Zeit
 * **/
template<int Di>
class MaxShiftFive {
public:
    enum {
        value = 5 * MaxShiftFive<Di - 1>::value
    };
};




template<>
class MaxShiftFive<1> {
public:
    enum {
        value = 5
    };
};


unsigned int MultiDimTwoCompass::maxShift =IntCases<DimensionSparseGrid>::value;

unsigned int MultiDimCompass::maxShift = MaxShift<DimensionSparseGrid>::value;




unsigned int MultiDimFiveCompass::maxShift = MaxShiftFive<DimensionSparseGrid>::value;



IndexDimension IndexDimension::nextTwoStep(MultiDimCompass *mc, Depth T) const {
// Wie normaler Kompass, aber nur zwei Schritte. Randnahe Punkte auslassen.

    IndexDimension J = (*this);
    while (mc->goon()) {

        for (int d = 0; d < DimensionSparseGrid; ++d) {
            int t = T.at(d);
            Richtung r = mc->getRichtung(d);

            if (r == Links) {
                if (J.nextLeft(d, t).isNotAtBoundary(d)) J = J.nextLeft(d, t).nextLeft(d, t);
                else { // break outer for loop
                    goto startagain;
                }
            }

            if (r == Rechts) {
                if (J.nextRight(d, t).isNotAtBoundary(d)) { J = J.nextRight(d, t).nextRight(d, t); }
                else { goto startagain; }
            }

        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        J = (*this);
        mc->operator++();
    }

    stop:

    return J;
}


IndexDimension IndexDimension::nextFive(MultiDimFiveCompass *mc, Depth T) const {

    IndexDimension J = (*this);
    IndexDimension Test;
    while (mc->goon()) {

        for (int d = 0; d < DimensionSparseGrid; ++d) {
            int t = T.at(d);
            Richtung r = mc->getRichtung(d);

            if (r == Links) {
                Test = J.nextLeft(d,t);
                if (Test.isNotAtBoundary(d)) J = Test;
                else // break outer for loop
                    goto startagain;
            }

            if (r == LinksLinks) {
                Test = J.nextLeft(d,t);
                if (Test.isNotAtBoundary(d)) J = Test.nextLeft(d, t);
                else goto startagain;
            }
            if (r == Rechts) {
                Test = J.nextRight(d,t);
                if (Test.isNotAtBoundary(d)) J = Test;
                else goto startagain;
            }
            if (r == RechtsRechts) {
                Test = J.nextRight(d,t);
                if (Test.isNotAtBoundary(d)) J = Test.nextRight(d, t);
                else goto startagain;
            }

        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        J = (*this);
        mc->operator++();
    }

    stop:

    return J;
}

void printDir(Richtung r) {
    if (r == Links) cout << "Links";
    if (r == LinksLinks) cout << "LinksLinks";
    if (r == Rechts) cout << "Rechts";
    if (r == RechtsRechts) cout << "RechtsRechts";
    if (r == Mitte) cout << "Mitte";
}

IndexDimension IndexDimension::nextFive(MultiDimFiveCompass *mc, Depth T, double *stencilvalue) const {
    double val = 1.0;
    IndexDimension J = (*this);
    IndexDimension TestIndex;

    while (mc->goon()) {
        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = T.at(d);
            Richtung r = mc->getRichtung(d);

            if (r == Links) {
                TestIndex = J.nextLeft(d,t);
                if (TestIndex.isNotAtBoundary(d)) {
                    val = -0.6 * val;
                    J = TestIndex;
                } else // break outer for loop
                    goto startagain;

            }

            if (r == LinksLinks) {
                TestIndex = J.nextLeft(d,t);
                if (TestIndex.isNotAtBoundary(d)) {
                    J = TestIndex.nextLeft(d, t);
                    val = 0.1 * val;
                } else goto startagain;
            }
            if (r == Rechts) {
                TestIndex = J.nextRight(d,t);
                if (TestIndex.isNotAtBoundary(d)) {
                    val = -0.6 * val;
                    J = TestIndex;
                } else goto startagain;
            }
            if (r == RechtsRechts) {
                TestIndex = J.nextRight(d,t);
                if (TestIndex.isNotAtBoundary(d)) {
                    J = TestIndex.nextRight(d, t);
                    val = 0.1 * val;
                } else goto startagain;
            }

            if (r == Mitte) {
                if (J.getDepth(d) > 1)
                    if (J.nextRight(d, t).isAtRightBoundary(d) || J.nextLeft(d, t).isAtLeftBoundary(d))
                        val = 0.9 * val;
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        val = 1.0;
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.
    if (mc->goon() == false && J == *this) val = 0.0;
    *stencilvalue = val;

    return J;
}


IndexDimension IndexDimension::nextFiveP(MultiDimFiveCompass *mc, Depth& T, double *stencilvalue) const {
   double val = 1.0;
    IndexDimension J = (*this);

   IndexDimension TestIndex;

   for (int d = 0; d < DimensionSparseGrid && mc->goon(); ++d) {

           int t = T.at(d);

           Richtung r = mc->getRichtung(d);

                     if (r == Links) {
                         TestIndex = J.nextLeft(d);
                         if (TestIndex.isNotAtBoundary(d)) {
                             val = -0.6 * val;
                             J = TestIndex;
                         } else {
                            mc->operator++();
                            d--;
                            continue;
                         }


                     }

                     if (r == LinksLinks) {
                         TestIndex = J.nextLeft(d);
                         if (TestIndex.isNotAtBoundary(d)) {
                             J = TestIndex.nextLeft(d, t);
                             val = 0.1 * val;
                         } else {
                             mc->operator++();
                             d--;
                             continue;
                         }
                     }
                     if (r == Rechts) {
                         TestIndex = J.nextRight(d);
                         if (TestIndex.isNotAtBoundary(d)) {
                             val = -0.6 * val;
                             J = TestIndex;
                         } else {
                             mc->operator++();
                             d--;
                             continue;
                         }
                     }
                     if (r == RechtsRechts) {
                         TestIndex = J.nextRight(d);
                         if (TestIndex.isNotAtBoundary(d)) {
                             J = TestIndex.nextRight(d, t);
                             val = 0.1 * val;
                         } else {
                             mc->operator++();
                             d--;
                             continue;
                         }
                     }

                     if (r == Mitte) {
                         if (J.getDepth(d) > 1)
                             if (J.isRandnah(d) || J.isRandnah(d))
                                 val = 0.9 * val;
                     }


       }


   //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.
   if (mc->goon() == false && J == *this) val = 0.0;
   *stencilvalue = val;

    return J;
}

IndexDimension IndexDimension::nextFive(MultiDimFiveCompass *mc, Depth T, double *stencilvalue, bool *last) const {
    double val = 1.0;
    IndexDimension J = (*this);

    while (mc->goon()) {
        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = T.at(d);
            Richtung r = mc->getRichtung(d);

            if (r == Links) {
                if (J.nextLeft(d, t).isNotAtBoundary(d)) {
                    val = -0.6 * val;
                    J = J.nextLeft(d, t);
                } else // break outer for loop
                    goto startagain;

            }

            if (r == LinksLinks) {
                if (J.nextLeft(d, t).isNotAtBoundary(d)) {
                    J = J.nextLeft(d, t).nextLeft(d, t);
                    val = 0.1 * val;
                } else goto startagain;
            }
            if (r == Rechts) {
                if (J.nextRight(d, t).isNotAtBoundary(d)) {
                    val = -0.6 * val;
                    J = J.nextRight(d, t);
                } else goto startagain;
            }
            if (r == RechtsRechts) {
                if (J.nextRight(d, t).isNotAtBoundary(d)) {
                    J = J.nextRight(d, t).nextRight(d, t);
                    val = 0.1 * val;
                } else goto startagain;
            }

            if (r == Mitte) {
                if (J.getDepth(d) > 1)
                    if (J.nextRight(d, t).isAtRightBoundary(d) || J.nextLeft(d, t).isAtLeftBoundary(d))
                        val = 0.9 * val;
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        val = 1.0;
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.
    if (mc->goon() == false && J == *this) {
        *last = false;
        val = 0.0;
    }
    *stencilvalue = val;

    return J;
}

IndexDimension IndexDimension::nextFive_Neumann(MultiDimFiveCompass *mc, Depth &T, double *stencilvalue) const {
    double val = 1.0;
    IndexDimension J = (*this);

    while (mc->goon()) {

        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = T.at(d);
            Richtung r = mc->getRichtung(d);
            //if(J.getDepth(d)==1) r = Mitte;

            if (r == Links) {
                if (J.getDepth(d) < 1) goto startagain;

                if (J.isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t);
                if (J.isAtLeftBoundary(d) && t > 1) val = -1.2 * val;
                else if (J.isAtLeftBoundary(d) && t == 1) val = -1.0 * val;
                else val = -0.6 * val;

            }

            if (r == LinksLinks) {
                if (J.getDepth(d) < 2) goto startagain;
                if (J.isAtLeftBoundary(d) || J.nextLeft(d, t).isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t).nextLeft(d, t);
                val = 0.1 * val;

            }

            if (r == Rechts) {
                if (J.getDepth(d) < 1) goto startagain;

                if (J.isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t);
                if (J.isAtRightBoundary(d) && t > 1) val = -1.2 * val;
                else if (J.isAtRightBoundary(d) && t == 1) val = -1.0 * val;
                else val = -0.6 * val;
            }
            if (r == RechtsRechts) {
                if (J.getDepth(d) < 2) goto startagain;

                if (J.isAtRightBoundary(d) || J.nextRight(d, t).isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t).nextRight(d, t);
                val = 0.1 * val;


            }

            if (r == Mitte) {


                if (J.getDepth(d) < 2) val = 1.0 * val;
                else {
                    if (J.nextLeft(d, t).isAtLeftBoundary(d) && !J.nextRight(d, t).isAtRightBoundary(d))
                        val = 1.1 * val;
                    if (J.nextRight(d, t).isAtRightBoundary(d) && !J.nextLeft(d, t).isAtLeftBoundary(d))
                        val = 1.1 * val;
                }
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        val = 1.0;
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.

    if (mc->goon() == false && J == *this) val = 0.0;

    *stencilvalue = val;

    return J;
}


IndexDimension IndexDimension::nextFive_Neumann(MultiDimFiveCompass *mc, Depth &T, double *stencilvalue, bool* last) const {
    double val = 1.0;
    IndexDimension J = (*this);

    while (mc->goon()) {

        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = T.at(d);
            Richtung r = mc->getRichtung(d);
            //if(J.getDepth(d)==1) r = Mitte;

            if (r == Links) {
                if (J.getDepth(d) < 1) goto startagain;

                if (J.isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t);
                if (J.isAtLeftBoundary(d) && t > 1) val = -1.2 * val;
                else if (J.isAtLeftBoundary(d) && t == 1) val = -1.0 * val;
                else val = -0.6 * val;

            }

            if (r == LinksLinks) {
                if (J.getDepth(d) < 2) goto startagain;
                if (J.isAtLeftBoundary(d) || J.nextLeft(d, t).isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t).nextLeft(d, t);
                val = 0.1 * val;

            }

            if (r == Rechts) {
                if (J.getDepth(d) < 1) goto startagain;

                if (J.isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t);
                if (J.isAtRightBoundary(d) && t > 1) val = -1.2 * val;
                else if (J.isAtRightBoundary(d) && t == 1) val = -1.0 * val;
                else val = -0.6 * val;
            }
            if (r == RechtsRechts) {
                if (J.getDepth(d) < 2) goto startagain;

                if (J.isAtRightBoundary(d) || J.nextRight(d, t).isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t).nextRight(d, t);
                val = 0.1 * val;


            }

            if (r == Mitte) {


                if (J.getDepth(d) < 2) val = 1.0 * val;
                else {
                    if (J.nextLeft(d, t).isAtLeftBoundary(d) && !J.nextRight(d, t).isAtRightBoundary(d))
                        val = 1.1 * val;
                    if (J.nextRight(d, t).isAtRightBoundary(d) && !J.nextLeft(d, t).isAtLeftBoundary(d))
                        val = 1.1 * val;
                }
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        val = 1.0;
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.

    if (mc->goon() == false && J == *this){
        *last = false;
        val = 0.0;
    }

    *stencilvalue = val;

    return J;
}


IndexDimension IndexDimension::nextThree_Stencil(MultiDimCompass *mc, Depth T, double *stencilvalue) const {
    double val = 1.0;
    IndexDimension J = (*this);

    while (mc->goon()) {

        for (int d = 0; d < DimensionSparseGrid; ++d) {

            int t = T.at(d);
            Richtung r = mc->getRichtung(d);
            //if(J.getDepth(d)==1) r = Mitte;



            if (r == Links) {
                if (J.getDepth(d) < 2) goto startagain;

                if (J.isAtLeftBoundary(d)) goto startagain;
                J = J.nextLeft(d, t);
                if (J.isAtLeftBoundary(d)) val = -1.2 * val;
                else val = -0.6 * val;

            }


            if (r == Rechts) {
                if (J.getDepth(d) < 2) goto startagain;

                if (J.isAtRightBoundary(d)) goto startagain;
                J = J.nextRight(d, t);
                if (J.isAtRightBoundary(d)) val = -1.2 * val;
                else val = -0.6 * val;
            }


            if (r == Mitte) {


                if (J.getDepth(d) < 2) val = 1.0 * val;
                else {
                    if (J.nextLeft(d, t).isAtLeftBoundary(d) && !J.nextRight(d, t).isAtRightBoundary(d))
                        val = 1.1 * val;
                    if (J.nextRight(d, t).isAtRightBoundary(d) && !J.nextLeft(d, t).isAtLeftBoundary(d))
                        val = 1.1 * val;
                }
            }


        }

        // forloop wurde für alle Dimensionen Erfolgreich ausgeführt -> gehe zu stop und return neuen Index
        goto stop;

        startagain:
        val = 1.0;
        J = (*this);
        mc->operator++();


    }

    stop:
    //Teste J darf nur == this sein, falls es über Mitte Mitte Mitte ... "generiert" wurde.

    if (mc->goon() == false && J == *this) val = 0.0;

    *stencilvalue = val;

    return J;
}


IndexDimension IndexDimension::Maximum(IndexDimension &A, IndexDimension &B) {
    IndexDimension P = A;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (P.coordinate(d) < B.coordinate(d)) P.replace(d, B.getIndex(d));
    }
    return P;
}

IndexDimension IndexDimension::Minimum(IndexDimension &A, IndexDimension &B) {
    IndexDimension P = A;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (P.coordinate(d) > B.coordinate(d)) P.replace(d, B.getIndex(d));
    }
    return P;
}


IndexDimension IndexDimension::Maximum(IndexDimension &A, IndexDimension &B, int *changes) {
    IndexDimension P = A;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (P.coordinate(d) < B.coordinate(d)) {
            P.replace(d, B.getIndex(d));
            changes[d] += 1;
        }
    }
    return P;

}

IndexDimension IndexDimension::Minimum(IndexDimension &A, IndexDimension &B, int *changes) {
    IndexDimension P = A;
    for (int d = 0; d < DimensionSparseGrid; ++d) {
        if (P.coordinate(d) > B.coordinate(d)) {
            P.replace(d, B.getIndex(d));
            changes[d] += 1;
        }
    }
    return P;
}
