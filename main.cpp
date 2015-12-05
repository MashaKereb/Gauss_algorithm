#include <fstream>
#include"gaussalgorithm.h"

using namespace std;

int main(){

    /*
    the intud data must be stored in the "input.txt" file
    in the folloving form:
        n m
        <matrixLine 1>  <freeMember 1>
        <matrixLine 2>  <freeMember 2>
        ................................
        <matrixLine m>  <freeMember m>

        where n = numberOfVariables, m = numberOfEquations
    */

    GaussAlgorithm ga;
    ifstream file("input.txt");
    ga.readFrom(file);

    int numberOfTreads = 0;
    cout << "Enter number of threads: ";
    cin >> numberOfTreads;

    while (numberOfTreads <= 0){
        cout << endl;
        cout << "This number must be nonegative! " << endl;
         cout << "Enter number of threads: ";
         cin >> numberOfTreads;

    }
    int numberOfSolutions = ga.findSolution(numberOfTreads);
    if (numberOfSolutions == 0) cout << "There is not any solution!";
    else {
        if(numberOfSolutions > 1) cout << "Thre are unlimit number of solutions!\n" <<
                                          "This is one of them: "<<endl;
    const double* solution = ga.getSolution();
    int n = ga.getNumberOfVariables();
        for(int i = 0; i < n; i++)
            cout << "x" << i << " = " << solution[i] << endl;

    }
    return 0;
}

