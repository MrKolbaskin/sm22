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
const int RECV_POINTS_TAG = 2;
const int SEND_POINTS_TAG = 3;

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
    srand(time(nullptr));
    int n = 0;
    double sumOfValues = 0, monteCarloValue = 0, tmpSum = 0;
    double points[3];

    while (abs(ANALITICAL_VALUE - monteCarloValue) > precision){
        points[0] = getRandomX();
        points[1] = getRandomY();
        points[2] = getRandomZ();
        
        MPI_Reduce(&tmpSum, MPI_IN_PLACE, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);

        n += numprocs;
        sumOfValues += tmpSum;

        monteCarloValue = VOLUME_P * (sumOfValues / double(numprocs));
    }

    cout << "Monte-Carlo Value: "<< monteCarloValue << endl;
    cout << "Diff Value: "<< abs(ANALITICAL_VALUE - monteCarloValue) << endl;
    cout << "Random points: "<< n << endl;
    cout << "Time: " << time(nullptr) - startTime << endl;
}

void workerLogic() {
    double points[3], value;
    int isActive = 1;

    while (isActive){
        MPI_Recv(&isActive, 1, MPI_INT, 0, ACTIVE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        if (isActive){
            MPI_Recv(points, 3, MPI_DOUBLE, 0, RECV_POINTS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            value = calcFunction(points[0], points[1], points[2]);

            MPI_Reduce(&value, nullptr, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
            // MPI_Send(value, 1, MPI_DOUBLE, 0, SEND_POINTS_TAG, MPI_COMM_WORLD);
        }
    }
}


int main(int argc, char *argv[]){
    time_t startTime = time(nullptr);
    
    if (argc < 2){
        cout << "You need specify precision" << endl;
        return 1;
    }

    double precision = atof(argv[1]);
    double reducedSum;

    int numprocs, myid;
    double tmpSum=0, sumOfValues=0;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    while (1) {
        
        MPI_Reduce(0, &tmpSum, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        if (myid){
            monteCarloValue = VOLUME_P * (sumOfValues / double(numprocs));

            if (abs(ANALITICAL_VALUE - monteCarloValue) > precision){
                break;
            }
        } else {
            MPI_Recv(&isActive, 1, MPI_INT, 0, ACTIVE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }


    MPI_Finalize();
    return 0;
}