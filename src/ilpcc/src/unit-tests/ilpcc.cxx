#include "andres/ilp/ilpcc.hxx"
#include "andres/ilp/cplex.hxx"

inline void test(const bool& x) { 
    if(!x) throw std::logic_error("test failed."); 
}

int main() {
    typedef andres::ilp::Cplex<double> Cplex;
    typedef andres::ilp::ILPCC<Cplex> ILPCC;

    ILPCC ilpcc;
    ilpcc.addVariables(6, ILPCC::Bool);

    // set objective function
    {
        double coefficients[] = {2, -2, -4, -3, 2, -2};
        ilpcc.setLinearCoefficients(coefficients);
    }

    ilpcc.optimize();
    test(ilpcc.value(0) == 0);
    test(ilpcc.value(1) == 1);
    test(ilpcc.value(2) == 1);
    test(ilpcc.value(3) == 1);
    test(ilpcc.value(4) == 0);
    test(ilpcc.value(5) == 1);

    // add constraint
    {
        size_t variableIndices[] = {0, 1, 2};
        double coefficients[] = {1, 1, 1};

        ilpcc.addConstraint(
            variableIndices, variableIndices + 3, 
            coefficients,
            0, 1
        );
    }

    ilpcc.optimize();
    test(ilpcc.value(0) == 0);
    test(ilpcc.value(1) == 0);
    test(ilpcc.value(2) == 1);
    test(ilpcc.value(3) == 1);
    test(ilpcc.value(4) == 0);
    test(ilpcc.value(5) == 1);

    // add constraint
    {
        size_t variableIndices[] = {3, 4, 5};
        double coefficients[] = {1, 1, 1};
        ilpcc.addConstraint(
            variableIndices, variableIndices + 3, 
            coefficients,
            0, 1
        );
    }

    ilpcc.optimize();
    test(ilpcc.value(0) == 0);
    test(ilpcc.value(1) == 0);
    test(ilpcc.value(2) == 1);
    test(ilpcc.value(3) == 1);
    test(ilpcc.value(4) == 0);
    test(ilpcc.value(5) == 0);

    // add constraint
    {
        size_t variableIndices[] = {2, 3};
        double coefficients[] = {1, 1};
        ilpcc.addConstraint(
            variableIndices, variableIndices + 2, 
            coefficients,
            0, 1
        );
    }

    ilpcc.optimize();
    test(ilpcc.value(0) == 0);
    test(ilpcc.value(1) == 0);
    test(ilpcc.value(2) == 1);
    test(ilpcc.value(3) == 0);
    test(ilpcc.value(4) == 0);
    test(ilpcc.value(5) == 1);

    return 0;
}
