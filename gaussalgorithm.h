#ifndef GAUSSALGORITHM_H
#define GAUSSALGORITHM_H

#include <thread>
#include <iostream>
#include <vector>
#include <stdexcept>

#define epsilon 0.0001
#define infinity 5000000

class GaussAlgorithm {
    int _numberOfVariables, _numberOfEquations;
    std::vector <double*> matrix;
    double* solution;

public:
    GaussAlgorithm();

    const double* getSolution () const;
    int getNumberOfVariables() const;


    void readFrom(std::istream& input);
    //returns the number of solutions
    int findSolution(int numberOfThreads = 1) throw(std::invalid_argument);
    void clear();

    ~GaussAlgorithm();
private:

    void iterate(std::vector<double*>, double* matrix, int currentPosition);
    void copyEquation(const double* from, double* to)const;





};

#endif // GAUSSALGORITHM_H
