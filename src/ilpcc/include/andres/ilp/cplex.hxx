#pragma once
#ifndef ANDRES_CPLEX_HXX
#define ANDRES_CPLEX_HXX

#include <stdexcept>
#include <sstream>
#include <vector>
#include <map>

#include "ilcplex/ilocplex.h"

namespace andres {
namespace ilp {

template<class T = double>
class Cplex {
public:
    typedef T Value;

    enum VariableType {Int, Float, Bool};

    enum PreSolver {
        PRE_SOLVER_AUTO, 
        PRE_SOLVER_PRIMAL, 
        PRE_SOLVER_DUAL, 
        PRE_SOLVER_NONE
    };

    enum LPSolver {
        LP_SOLVER_AUTO, 
        LP_SOLVER_PRIMAL_SIMPLEX, 
        LP_SOLVER_DUAL_SIMPLEX, 
        LP_SOLVER_NETWORK_SIMPLEX, 
        LP_SOLVER_BARRIER, 
        LP_SOLVER_SIFTING, 
        LP_SOLVER_CONCURRENT
    };

    Cplex();
    ~Cplex();
    void operator=(const Cplex<T>&);
    void clear();

    // access properties
    size_t numberOfVariables() const; 
    Value value(const size_t) const;
    size_t numberOfThreads() const;
    Value absoluteGap() const;
    Value relativeGap() const;
    size_t numberOfConstraints() const;

    // manipulate properties
    void setNumberOfThreads(const size_t);
    void setAbsoluteGap(const Value);
    void setRelativeGap(const Value);
    void setVerbosity(const bool);
    void setLPSolver(const LPSolver);
    void setPreSolver(const PreSolver, const int = -1);

    // set up a problem
    size_t addVariable(const VariableType); 
    size_t addVariables(const size_t, const VariableType);
    void setVariableType(const size_t, const VariableType);
    void setLinearCoefficient(const size_t, const Value);
    template<class CoefficientIterator>
        void setLinearCoefficients(CoefficientIterator);
    template<class VariableIterator, class CoefficientIterator>
        void setLinearCoefficients(VariableIterator, VariableIterator,
            CoefficientIterator);
    void setQuadraticCoefficient(const size_t, const size_t, const Value);
    void addConstraint(const size_t, const Value, const Value); 
    template<class VariableIndexIterator, class CoefficientIterator>
        void addConstraint(VariableIndexIterator, VariableIndexIterator,
                           CoefficientIterator, const Value, const Value);

    // TODO: write function setObjectiveFunction(ObjectiveFunction)
    void addLinearCoefficient(const size_t, const Value);
    void addQuadraticCoefficient(const size_t, const size_t, const Value);
    void setObjectiveFunction();

    // solve a problem
    template<class Iterator>
        void setStart(Iterator);
    void registerStart(const size_t, const double);
    void setRegisteredStart();
    Value optimize();
  
private:
    IloNumVar::Type cplexVariableType(const VariableType) const; 

    IloEnv ilpEnvironment_;
    IloObjective ilpObjective_;
    IloExpr ilpExpression_;
    IloModel ilpModel_;
    IloCplex ilpSolver_;
    IloNumArray ilpSolution_;
    IloNumVarArray ilpVariables_;
    IloNumArray ilpStartValues_;
    std::map<size_t, IloConversion> ilpVariableConversions_;
    IloRangeArray constraints_; 
};

template<class T>
inline
Cplex<T>::Cplex()
:   ilpEnvironment_(),
    ilpObjective_(ilpEnvironment_),
    ilpExpression_(ilpEnvironment_),
    ilpModel_(ilpEnvironment_),
    ilpSolver_(ilpModel_),
    ilpSolution_(ilpEnvironment_),
    ilpVariables_(ilpEnvironment_),
    ilpStartValues_(ilpEnvironment_),
    ilpVariableConversions_(),
    constraints_(ilpEnvironment_)
{
    setAbsoluteGap(0);
    setRelativeGap(0);
    setVerbosity(false);
    try {
        ilpObjective_.setSense(IloObjective::Minimize);
        ilpModel_.add(ilpObjective_);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline
Cplex<T>::~Cplex() {
    ilpEnvironment_.end(); // release memory
}

// TODO: test this operator
template<class T>
inline void
Cplex<T>::operator=(
    const Cplex<T>& other
) {
    throw std::runtime_error("Internal error. Assignment operator not implemented.");
}

template<class T>
inline void
Cplex<T>::clear() {
    ilpSolver_.clear();
    ilpSolver_.clearModel();
    // ilpModel_.???
    // ilpObjective_.???    
    ilpExpression_.clear();
    ilpSolution_.clear();
    ilpVariables_.clear();
    ilpStartValues_.clear();
    ilpVariableConversions_.clear();
    constraints_.clear();

    try {
        ilpSolver_.extract(ilpModel_);
        ilpObjective_.setSense(IloObjective::Minimize);
        ilpModel_.add(ilpObjective_);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline void
Cplex<T>::setNumberOfThreads (
    const size_t numberOfThreads
) {
    try {
        ilpSolver_.setParam(IloCplex::Threads, numberOfThreads);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline void
Cplex<T>::setAbsoluteGap(
    const T gap
) {
    try {
        ilpSolver_.setParam(IloCplex::EpAGap, gap);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline void
Cplex<T>::setRelativeGap(
    const T gap
) {
    try {
        ilpSolver_.setParam(IloCplex::EpGap, gap);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline void
Cplex<T>::setVerbosity(
    const bool verbosity
) {
    try {
        if(verbosity) {
            ilpSolver_.setParam(IloCplex::MIPDisplay, 2);
            ilpSolver_.setParam(IloCplex::BarDisplay, 1);
            ilpSolver_.setParam(IloCplex::NetDisplay, 1);
            ilpSolver_.setParam(IloCplex::SiftDisplay, 1);
            ilpSolver_.setParam(IloCplex::SimDisplay, 1);
        }
        else {
            ilpSolver_.setParam(IloCplex::MIPDisplay, 0);
            ilpSolver_.setParam(IloCplex::BarDisplay, 0);
            ilpSolver_.setParam(IloCplex::NetDisplay, 0);
            ilpSolver_.setParam(IloCplex::SiftDisplay, 0);
            ilpSolver_.setParam(IloCplex::SimDisplay, 0);
        }
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline void
Cplex<T>::setPreSolver(
    const PreSolver preSolver,
    const int passes 
) {
    try {
        switch(preSolver) {
        case PRE_SOLVER_NONE:
            ilpSolver_.setParam(IloCplex::PreInd, CPX_OFF);
            ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_OFF);
            return;
        case PRE_SOLVER_AUTO:
            ilpSolver_.setParam(IloCplex::PreInd, CPX_ON);
            ilpSolver_.setParam(IloCplex::PreDual, 0); // default
            ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_ON);
            break;
        case PRE_SOLVER_PRIMAL:
            ilpSolver_.setParam(IloCplex::PreInd, CPX_ON);
            ilpSolver_.setParam(IloCplex::PreDual, -1); 
            ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_ON);
            break;
        case PRE_SOLVER_DUAL:
            ilpSolver_.setParam(IloCplex::PreInd, CPX_ON);
            ilpSolver_.setParam(IloCplex::PreDual, 1); 
            ilpSolver_.setParam(IloCplex::RelaxPreInd, CPX_ON);
            break;
        }
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline void
Cplex<T>::setLPSolver(
    const LPSolver lpSolver
) {
    try {
        switch(lpSolver) {
        case LP_SOLVER_AUTO:
            ilpSolver_.setParam(IloCplex::RootAlg, 0); // default
            ilpSolver_.setParam(IloCplex::NodeAlg, 0); // default
            break;
        case LP_SOLVER_PRIMAL_SIMPLEX:
            ilpSolver_.setParam(IloCplex::RootAlg, 1);
            ilpSolver_.setParam(IloCplex::NodeAlg, 1);
            break;
        case LP_SOLVER_DUAL_SIMPLEX:
            ilpSolver_.setParam(IloCplex::RootAlg, 2);
            ilpSolver_.setParam(IloCplex::NodeAlg, 2);
            break;
        case LP_SOLVER_NETWORK_SIMPLEX:
            ilpSolver_.setParam(IloCplex::RootAlg, 3);
            ilpSolver_.setParam(IloCplex::NodeAlg, 3);
            break;
        case LP_SOLVER_BARRIER:
            ilpSolver_.setParam(IloCplex::RootAlg, 4);
            ilpSolver_.setParam(IloCplex::NodeAlg, 4);
            break;
        case LP_SOLVER_SIFTING:
            ilpSolver_.setParam(IloCplex::RootAlg, 5);
            ilpSolver_.setParam(IloCplex::NodeAlg, 5);
            break;
        case LP_SOLVER_CONCURRENT:
            ilpSolver_.setParam(IloCplex::RootAlg, 6);
            ilpSolver_.setParam(IloCplex::NodeAlg, 6);
            break;
        }
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Add one variable to an ILP (delayed column generation).
///
/// \return Index of the newly inserted variable.
///
template<class T>
inline size_t
Cplex<T>::addVariable(
    const VariableType variableType
) {
    return addVariables(1, variableType);
}

/// Add variables to an ILP (delayed column generation).
///
/// \param number Number of variables.
/// \param cit Iterator pointing to the beginning of a sequence of coefficients.
/// \return Index of the first newly inserted variable.
///
template<class T>
inline size_t
Cplex<T>::addVariables(
    const size_t number,
    const VariableType variableType
) {
    try {
        size_t offset = ilpVariables_.getSize();
        ilpVariables_.add(
            IloNumVarArray(
                ilpEnvironment_, 
                number, 
                -IloInfinity, IloInfinity, 
                cplexVariableType(variableType)
            )
        );
        ilpStartValues_.add(IloNumArray(ilpEnvironment_, number));
        for(size_t j = 0; j < number; ++j) {
            ilpStartValues_[j] = 1;			// test
        }
        return offset;
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Change the type of a variable.
///
/// variable Index of the variable.
/// variableType VariableType.
///
/// This is a 'no-worries implementation' of the conversion mechanism of
/// IBM ILOG Cplex. A conversions will be added or removed as necessary.
///
template<class T>
inline void
Cplex<T>::setVariableType(
    const size_t variable, 
    const VariableType variableType
) {
    try {
        std::map<size_t, IloConversion>::iterator it 
            = ilpVariableConversions_.find(variable);
        if(ilpVariables_[variable].getType() == cplexVariableType(variableType)) { 
        // if the target type is the original type of the variable
            if(it != ilpVariableConversions_.end()) { 
            // if a conversion of this variable exists
                // remove conversion
                it->second.end();
                ilpVariableConversions_.erase(it);
            }
        }
        else { // if the target type is not the original type of the variable
            if(it == ilpVariableConversions_.end()) {
            // if no conversion of the variable exists
                // insert new conversion
                ilpVariableConversions_[variable] = IloConversion(
                    ilpEnvironment_, 
                    ilpVariables_[variable],
                    cplexVariableType(variableType)
                );
                ilpModel_.add(ilpVariableConversions_[variable]);
            }
            else {
            // if a conversion of the variable exists
                // remove conversion
                it->second.end();
                ilpVariableConversions_.erase(it);
                // insert new conversion
                ilpVariableConversions_[variable] = IloConversion(
                    ilpEnvironment_, 
                    ilpVariables_[variable],
                    cplexVariableType(variableType)
                );
                ilpModel_.add(ilpVariableConversions_[variable]);
            }
        }
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Set the linear coefficient of one variable.
///
/// WARNING: Calling this function repeatedly, for a large number of varibales,
/// can be slow. Consider using setLinearCoefficients() instead.
///
/// \sa setLinearCoefficients()
///
template<class T>
inline void
Cplex<T>::setLinearCoefficient(
    const size_t variable,
    const Value coefficient
) {
    try {
        ilpObjective_.setLinearCoef(ilpVariables_[variable], coefficient);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Set the linear coefficients of all variables.
///
template<class T>
template<class CoefficientIterator>
inline void
Cplex<T>::setLinearCoefficients(
    CoefficientIterator cit
) {
    try {
        IloNumArray coefficients(ilpEnvironment_, numberOfVariables());
        for(size_t j = 0; j < numberOfVariables(); ++j, ++cit) {
			coefficients[j] = *cit;
        }
        ilpObjective_.setLinearCoefs(ilpVariables_, coefficients);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Set the linear coefficients of a subset of variables.
///
template<class T>
template<class VariableIterator, class CoefficientIterator>
inline void
Cplex<T>::setLinearCoefficients(
    VariableIterator vit,
    VariableIterator vend,
    CoefficientIterator cit
) {
    for(; vit != vend; ++vit, ++cit) {
        setLinearCoefficient(*vit, *cit);
    }
}

// TODO: write doc
template<class T>
inline void
Cplex<T>::setQuadraticCoefficient(
    const size_t variable0,
    const size_t variable1,
    const Value coefficient
) {    
    try {
        ilpObjective_.setQuadCoef(
            ilpVariables_[variable0], 
            ilpVariables_[variable1], 
            coefficient
        );
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

// TODO: write doc
template<class T>
inline void
Cplex<T>::addLinearCoefficient(
    const size_t variable,
    const Value coefficient
) {    
    try {
        ilpExpression_ += coefficient * ilpVariables_[variable];
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

// TODO: write doc
template<class T>
inline void
Cplex<T>::addQuadraticCoefficient(
    const size_t variable0,
    const size_t variable1,
    const Value coefficient
) {   
    try {
        ilpExpression_ += coefficient * ilpVariables_[variable0] * ilpVariables_[variable1];
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

// TODO: write doc
template<class T>
inline void
Cplex<T>::setObjectiveFunction() {
    try {
        ilpObjective_.setExpr(ilpExpression_);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Return the number of variables.
///
template<class T>
inline size_t
Cplex<T>::numberOfVariables() const {
    try {
        return ilpVariables_.getSize();
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Return the value of a variable.
///
template<class T>
inline T
Cplex<T>::value (
    const size_t variableIndex
) const {
    try {
        return ilpSolution_[variableIndex];
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}


/// Return the number of threads used by IBM ILOG Cplex.
///
template<class T>
inline size_t
Cplex<T>::numberOfThreads() const {
    try {
        return ilpSolver_.getParam(IloCplex::Threads);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Return the absolute gap used as a stopping criterion.
///
template<class T>
inline T
Cplex<T>::absoluteGap() const {
    try {
        return ilpSolver_.getParam(IloCplex::EpAGap);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Return the relative gap used as a stopping criterion.
///
template<class T>
inline T
Cplex<T>::relativeGap() const {
    try {
        return ilpSolver_.getParam(IloCplex::EpGap);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Get the number of constraints.
///
template<class T>
inline size_t
Cplex<T>::numberOfConstraints() const {
	return constraints_.getSize();
}

/// Add a constraint.
///
/// \param vit Iterator pointing to the beginning of a sequence of variable indices.
/// \param vend Iterator pointing to the end of a sequence of variable indices.
/// \param coefficient Iterator pointing to the beginning of a sequence of coefficients
/// \param lowerBound Lower bound on the linear combination.
/// \param upperBound Upper bound on the linear combination.
///
template<class T>
template<class VariableIndexIterator, class CoefficientIterator>
inline void 
Cplex<T>::addConstraint(
    VariableIndexIterator vit,
    VariableIndexIterator vend,
    CoefficientIterator coefficient,
    const Value lowerBound,
    const Value upperBound
) {
    try {
        /*
        IloRange constraint(ilpEnvironment_, lowerBound, upperBound);
        for(; vit != vend; ++vit, ++coefficient) {
           constraint.setLinearCoef(ilpVariables_[*vit], *coefficient);             
        }
        ilpModel_.add(constraint);
        */

        constraints_.add(IloRange(ilpEnvironment_, lowerBound, upperBound));
        IloRange& constraint = constraints_[constraints_.getSize() - 1];
        for(; vit != vend; ++vit, ++coefficient) {
           constraint.setLinearCoef(ilpVariables_[*vit], *coefficient);             
        }
        ilpModel_.add(constraint);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Add a constraint on a single variable.
///
/// \param variableIndex Index of the variable.
/// \param lowerBound Lower bound on the value of the variable.
/// \param upperBound Upper bound on the value of the variable.
///
template<class T>
inline void 
Cplex<T>::addConstraint(
    const size_t variableIndex, 
    const Value lowerBound, 
    const Value upperBound
) {
    try {
        /*
        IloRange constraint(ilpEnvironment_, lowerBound, upperBound);
        constraint.setLinearCoef(ilpVariables_[variableIndex], 1);
        ilpModel_.add(constraint);
        */

        constraints_.add(IloRange(ilpEnvironment_, lowerBound, upperBound));
        IloRange& constraint = constraints_[constraints_.getSize() - 1];
        constraint.setLinearCoef(ilpVariables_[variableIndex], 1);
        ilpModel_.add(constraint);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

/// Set a start for the solver and delete existing ones (if any).
///
template<class T>
template<class Iterator>
inline void
Cplex<T>::setStart(
    Iterator valueIterator
) {
    try {
        // delete existing mip starts
        if(ilpSolver_.getNMIPStarts() != 0) {
            ilpSolver_.deleteMIPStarts(0, ilpSolver_.getNMIPStarts());
        }
        // add new mip start
        for(size_t j = 0; j < static_cast<size_t>(ilpVariables_.getSize()); ++j, ++valueIterator) {
            ilpStartValues_[j] = *valueIterator;
        }
        ilpSolver_.addMIPStart(ilpVariables_, ilpStartValues_);
    }
    catch(const IloException& err) {
        throw std::runtime_error(err.getMessage());
    }
}

// TODO: maybe clean up (not all variable may be set)
template<class T>
inline void
Cplex<T>::registerStart(
    const size_t variableIndex,
    const double value
) {
    ilpStartValues_[variableIndex] = value;
}

// TODO: maybe clean up (not all variable may be set)
template<class T>
inline void
Cplex<T>::setRegisteredStart() {
    ilpSolver_.addMIPStart(ilpVariables_, ilpStartValues_);
}

/// Run the IBM ILOG Cplex solver.
///
template<class T>
inline typename Cplex<T>::Value
Cplex<T>::optimize() {
    try {
        ilpSolver_.solve();
        ilpSolver_.getValues(ilpSolution_, ilpVariables_);
        return static_cast<Value>(ilpSolver_.getObjValue());
    }
    catch(const IloException& err) {
		throw std::runtime_error(err.getMessage());
    }
}

template<class T>
inline IloNumVar::Type
Cplex<T>::cplexVariableType(
    const VariableType variableType
) const {
    switch(variableType) {
    case Float:
        return IloNumVar::Float;
    case Int:
        return IloNumVar::Int;
    default: // Bool
        return IloNumVar::Bool;
    }
}

} // namespace ilp
} // namespace andres

#endif // #ifndef ANDRES_CPLEX_HXX
