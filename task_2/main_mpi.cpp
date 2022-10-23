#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include "mpi.h"

using namespace std;

const double MIN_X = -1;
const double MAX_X = 1;
const double MIN_Y = -1;
const double MAX_Y = 1;

const double ANALITICAL_VALUE = M_PI / 6.0;
const double VOLUME_P = 4;

const int ROOT = 0;

const int ACTIVE_TAG = 1;
const int POINTS_TAG = 2;

double getRandomNumber(){
    return rand() / double(RAND_MAX); 
}

double getRandomZ(){
    return getRandomNumber();
}

double getRandomX(){
    return MIN_X + (MAX_X - MIN_X) * getRandomNumber();
}

double getRandomY(){
    return MIN_Y + (MAX_Y - MIN_Y) * getRandomNumber();
}

double calcFunction(const double &x, const double &y, const double &z){
    if (x * x + y * y <= z * z) return sqrt(x * x + y * y);

    return 0;
}

void masterLogic(const double precision, const int numprocs) {
    time_t startTime = time(nullptr);
    srand(23);
    int n = 0, isActive = 1;
    double sumOfValues = 0, monteCarloValue = 0, tmpSum = 0, innerSum = 0;
    double points[3];

    while (abs(ANALITICAL_VALUE - monteCarloValue) > precision){
        MPI_Bcast(&isActive, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

        for (int i = 1; i < numprocs; i++){
            points[0] = getRandomX();
            points[1] = getRandomY();
            points[2] = getRandomZ();

            MPI_Send(points, 3, MPI_DOUBLE, i, POINTS_TAG, MPI_COMM_WORLD);
        }

        MPI_Reduce(&innerSum, &tmpSum, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        n += numprocs - 1;
        sumOfValues += tmpSum;

        monteCarloValue = VOLUME_P * (sumOfValues / double(n));
    }

    isActive = 0;
    MPI_Bcast(&isActive, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    cout << "Monte-Carlo Value: "<< monteCarloValue << endl;
    cout << "Diff Value: "<< abs(ANALITICAL_VALUE - monteCarloValue) << endl;
    cout << "Random points: "<< n << endl;
}

void workerLogic() {
    double points[3], value, outValue;
    int isActive = 1;

    while (isActive){
        MPI_Bcast(&isActive, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        
        if (isActive){
            MPI_Recv(points, 3, MPI_DOUBLE, 0, POINTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            value = calcFunction(points[0], points[1], points[2]);

            MPI_Reduce(&value, &outValue, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        }
    }
}


int main(int argc, char *argv[]){
    if (argc < 2){
        cout << "You need specify precision" << endl;
        return 1;
    }
    
    MPI_Init(&argc, &argv);
    int numprocs, myid;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    double startTime = MPI_Wtime();

    double precision = atof(argv[1]);
    double tmpSum=0, sumOfValues=0, reducedSum;

    if (!myid){
        masterLogic(precision, numprocs);
    } else {
        workerLogic();
    }

    double myTotalTime = MPI_Wtime() - startTime, totalTime;
    MPI_Reduce(&myTotalTime, &totalTime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
    if (!myid){
        cout << "Time: " << totalTime << endl;
    }

    MPI_Finalize();
    return 0;
}