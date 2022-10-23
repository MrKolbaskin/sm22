#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;

const double MIN_X = -1;
const double MAX_X = 1;
const double MIN_Y = -1;
const double MAX_Y = 1;

const double ANALITICAL_VALUE = M_PI / 6.0;
const double VOLUME_P = 4;

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

double calcFunction(const double x, const double y, const double z){
    if (x * x + y * y <= z * z) return sqrt(x * x + y * y);

    return 0;
}


int main(int argc, char *argv[]){
    time_t startTime = time(nullptr);
    
    if (argc < 2){
        cout << "You need specify precision" << endl;
        return 1;
    }

    double precision = atof(argv[1]);

    srand(23);
    int n = 0;
    double sumOfValues = 0;
    double monteCarloValue = 0;

    while (abs(ANALITICAL_VALUE - monteCarloValue) > precision){
        double x = getRandomX();
        double y = getRandomY();
        double z = getRandomZ();

        sumOfValues += calcFunction(x, y, z);
        n++;

        monteCarloValue = VOLUME_P * (sumOfValues / double(n));
    }


    cout << "Monte-Carlo Value: "<< monteCarloValue << endl;
    cout << "Diff Value: "<< abs(ANALITICAL_VALUE - monteCarloValue) << endl;
    cout << "Random points: "<< n << endl;
    cout << "Time: " << time(nullptr) - startTime << endl;

    return 0;
}