#include "functions.cpp"
using namespace std;



int main() {

    int nProcs, iProc;
    bool did_work;

    // SET THE BOUNDS FOR THE DATA
    int begin = -1000;
    int finish = 1000;

    // SET THE SIZE OF grid 1
    int start1 = -5;
    int end1 = 7;
    double interval1 = 0.5;
    int size1 = abs(end1 - start1);

    // SET THE SIZE OF grid 2
    int start2 = -2;
    int end2 = 6;
    double interval2 = 1;
    int size2 = abs(end2 - start2);

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &iProc);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

    // randomly calculate a seed for every processor
    double seed = random() % (((iProc + 1) * 1000000) - (iProc * 1000000) + 1) + (iProc * 1000000);
    srand(seed);

    // calculates what processor takes care of how much of the array
    int subSize1 = size1 / nProcs;
    if (size1 % nProcs != 0) ++subSize1;

    int subSize2 = size2 / nProcs;
    if (size2 % nProcs != 0) ++subSize2;


    // generates grid 1 (each processor handles a part of grid 1)
    map<double, double> grid1 = genArray(iProc, subSize1, size1, interval1, start1, end1, begin, finish);


    // finds what processor handles what part of the grid and stores it in array
    // prints array
    MPI_Barrier(MPI_COMM_WORLD);
    map<double, double> array1 = findLocations(&grid1, size1, interval1, start1, iProc);
    if (iProc == 0) {
        cout << "\n-------------------- ARRAYS --------------------\n\n" << "ARRAY ONE" << endl;
        map<double, double>::iterator it = array1.begin();
        while (it != array1.end()) {
            cout << "x = " << it->first << ", value = " << it->second << endl;
            ++it;
        }
    }
    
    

    // generates grid 2 (each processor handles a part of grid 2)
    map<double, double> grid2 = genArray(iProc, subSize2, size2, interval2, start2, end2, begin, finish);

    // prints array for second grid 
    MPI_Barrier(MPI_COMM_WORLD);
    map<double, double> array2 = findLocations(&grid2, size2, interval2, start2, iProc);
    MPI_Barrier(MPI_COMM_WORLD);
    if (iProc == 0) {
        cout << "\nARRAY TWO" << endl;
        map<double, double>::iterator it = array2.begin();
        while (it != array2.end()) {
            cout << "x = " << it->first << ", value = " << it->second << endl;
            ++it;
        }
    }


    if (iProc == 0) cout << "\n-------------------- GRIDS --------------------\n" << endl;

    // prints out both grids 1 and 2
    MPI_Barrier(MPI_COMM_WORLD);

    
    
    if (iProc == 0) cout << "GRID ONE" << endl;
    sleep(0.5);
    printGrid(&array1, &grid1, iProc);

    

    sleep(0.5);
    if (iProc == 0) cout << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    if (iProc == 0) cout << "GRID TWO" << endl;
    sleep(0.5);
    printGrid(&array2, &grid2, iProc);



    sleep(0.75);
    map<double, double> coeff1 = findCoeff(&array1, &array2);
    map<double, double> coeff2 = findCoeff(&array2, &array1);
    if (iProc == 0) {
        cout << "\n-------------------- testing coefficient --------------------\n" << endl;
        cout << "COEFFICIENT ONE" << endl;
        map<double, double>::iterator it = coeff1.begin();
        while (it != coeff1.end()) {
            cout << "x = " << it->first << ", coeff = " << it->second << endl;
            ++it;
        }

        cout << "\nCOEFFICIENT TWO" << endl;
        it = coeff2.begin();
        while (it != coeff2.end()) {
            cout << "x = " << it->first << ", coeff = " << it->second << endl;
            ++it;
        }
    }



    sleep(0.75);
    if (iProc == 0) {
        cout << "\n-------------------- testing message passing --------------------\n" << endl;
        cout << "TEST getValue()" << endl;
    }


    MPI_Barrier(MPI_COMM_WORLD);
    sleep(0.5);


    // CHANGE THIS VALUE TO TEST DIFFERENT X-POSITIONS
    double xPos = 0;
    // CHANGE THIS VALUE TO TEST DIFFERENT X-POSITIONS


    // testing the message passing interface
    double test;
    int procReceive = findProc(&array1, xPos);
    int procSend = findProc(&array2, xPos);
    getValue(&array1, &array2, &grid2, &coeff1, iProc, xPos, &test);
    if (iProc == procReceive) {
        cout << "x-position = " << xPos << endl;
        cout << "receiving processor = " << procReceive << endl;
        cout << "sending processor = " << procSend << endl;
        cout << "actual = " << test << endl;
    }
    sleep(0.5);
    if (iProc == procSend) cout << "expected = " << grid2.at(xPos) << endl;


    MPI_Barrier(MPI_COMM_WORLD);
    getValue(&array2, &array1, &grid1, &coeff2, iProc, xPos, &test);
    if (iProc == procSend) {
        cout << "\nTEST 2 getValue()" << endl;
        cout << "x-position = " << xPos << endl;
        cout << "receiving processor = " << procReceive << endl;
        cout << "sending processor = " << procSend << endl;
        cout << "actual = " << test << endl;
    }
    sleep(0.5);
    if (iProc == procReceive) cout << "expected = " << grid1.at(xPos) << endl;




    xPos = 1.5;
    test = 0;

    MPI_Barrier(MPI_COMM_WORLD);
    sleep(0.5);
    procReceive = getValue(&array1, &array2, &grid2, &coeff1, iProc, xPos, &test);
    if (iProc == procReceive) {
        cout << "\nTEST INTERPOLATION getValue()" << endl;
        cout << "x-position = " << xPos << endl;
        cout << "receiving processor = " << procReceive << endl;
        cout << "actual = " << test << endl;
    }



    MPI_Finalize();
    return 0;
}