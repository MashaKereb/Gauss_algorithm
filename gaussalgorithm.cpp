
#include "gaussalgorithm.h"
using namespace std;

GaussAlgorithm::GaussAlgorithm(){
    _numberOfEquations = 0;
    _numberOfVariables = 0;
    solution = nullptr;
}

void GaussAlgorithm::readFrom(istream& in){

    in >> _numberOfVariables;
    in >> _numberOfEquations;

    double* currEquation;
    for (int i = 0; i < _numberOfEquations; i++){

        currEquation = new double[_numberOfVariables + 1];

        for(int j = 0; j <= _numberOfVariables; j++)
            in >> currEquation[j];

        matrix.push_back(currEquation);
    }
}

const double* GaussAlgorithm::getSolution () const{
    return solution;
}

int GaussAlgorithm::getNumberOfVariables() const{
    return _numberOfVariables;
}

void GaussAlgorithm::iterate(vector<double*> matrix, double* currEquation, int position){

    int size = matrix.size();
    double multiplier = 0.0;

    for(int i = 0; i < size; i++){
         multiplier = matrix[i][position] / currEquation[position];

        for(int j = 0; j <= _numberOfVariables; j++){

            matrix[i][j] -= currEquation[j] * multiplier;
        }
    }

}

void GaussAlgorithm::copyEquation(const double* from, double* to)const{

    for (int i = 0; i <= _numberOfVariables; i++)
        to[i] = from[i];
}


int GaussAlgorithm::findSolution(int numberOfTreads) throw(invalid_argument) {
    if(numberOfTreads <= 0) throw invalid_argument("Incorect number of threads");
    if(numberOfTreads > _numberOfEquations) numberOfTreads = _numberOfEquations;



    thread threads[numberOfTreads];
    vector<double*> Equations[numberOfTreads];

    /* Matrix division into small groups beetween  the threads.
       Every next equation is cyclicaly added to the other group
    */
    for(int i = 0; i < _numberOfEquations; i++){
        Equations[i % numberOfTreads].push_back(matrix[i]);
    }

    /*
        in "positions[i]" will be stored information
        about value of which variable can be found from i-th
        equation in matrix.

        In the best way (if matrix are square, her lines are linearly independent)
        we will get diagonal matrix after iterations, and this information will be redandant.
        But when something goes wrong (current line is empty or current element is equal to zero)
        we need to change order of lines.
        But in this situation it is easier to change current variable, not whole line
        (it is also easier to check whether the system  has solution at all).
        So in this situation we need to save the order of variables.
    */
    vector<int>positions (_numberOfEquations);
    for(int i = 0; i < _numberOfEquations; i++)
        positions[i] = i;



    /*
        Here is the Gauss' algorithm!
        For each iteration every line in the matrix is taken,
        checked for being empty and passed
        into each thread for further processing.
    */



    double* currentEquation = new double[_numberOfVariables + 1]; //copy of current line of matrix
    double curr = 0.0;                                            //copy of current element
                            /* we need a copy, becouse information may be lost during iterations
                            */

    // these elements are used for checking
    // whether the system  has solution at all
    bool isEmpty = false;
    bool hasSolution = true;
    int _numberOfEmptyLines = 0;


    for(int i = 0; i < _numberOfEquations; i++){

        copyEquation( matrix[i], currentEquation);

        //here is the current line's checking
        //current element cannot be equal to zero,
        //so we need to find any other nonzero element
        //if it is possible
        if(::abs(currentEquation[positions[i]]) < epsilon) {
            int j = i + 1;
            while(::abs(currentEquation[positions[j]]) < epsilon){
                if(j == _numberOfVariables)
                    break;
                j++;
            }
            if(j == _numberOfVariables){
                if(::abs(currentEquation[j]) < epsilon) isEmpty = true;
                else hasSolution = false;
            }
            else std::swap(positions[i], positions[j]);
        }
        //if system has no solution we can and the algorithm
        if(hasSolution == false) {
            delete currentEquation;
            return 0;
        }
        //if line is empty (all elements are equal to zero)
        //we can miss it. But we need to mark that none of variables can be
        //fond from this equation
        if(isEmpty){
            isEmpty = false;
            _numberOfEmptyLines++;
            for(int j = _numberOfEquations - 1;j > i; j--)
                positions[j] = positions[j - 1];
            positions[i] -1;
            continue;
        }



        for(int j = 0; j < numberOfTreads; j++)
         threads[j] = thread(&GaussAlgorithm::iterate, this, Equations[j], currentEquation, positions[i]);

        for(int j = 0; j < numberOfTreads; j++)
         threads[j].join();

        //separately convert the current line to the appropriate form
        curr = currentEquation[positions[i]];
        for(int j = 0; j <=_numberOfVariables; j++){
         currentEquation[j] /= curr;
       }

        copyEquation(currentEquation, matrix[i]);
    }

    //filling solution's vector
     solution = new double[_numberOfVariables];
     for(int i = 0; i < _numberOfVariables; i++)
         solution[i] = 0.0;

    for(int i = 0; i < _numberOfEquations; i++) {
        if(positions[i] == -1) continue;
        solution[positions[i]] = matrix[i][_numberOfVariables];
    }
    //function returns the number of solutions
     if(_numberOfEquations - _numberOfEmptyLines == _numberOfVariables)
         return 1;
     return infinity;
}

void GaussAlgorithm::clear(){
   int size = matrix.size();
   for(int i = 0; i < size; i++)
       delete matrix[i];
   matrix.clear();
   delete solution;
   solution = nullptr;
   _numberOfEquations = 0;
   _numberOfVariables = 0;
}
GaussAlgorithm::~GaussAlgorithm(){
   int size = matrix.size();
   for(int i = 0; i < size; i++)
       delete matrix[i];
   delete solution;
}
