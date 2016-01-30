#pragma once
#ifndef ANDRES_ILP_ILPCC_HXX
#define ANDRES_ILP_ILPCC_HXX

#include <cassert>
#include <stdexcept>
#include <limits>
#include <vector>
#include <set>
#include <algorithm> // std::distance
#include <iostream>

#include "andres/partition.hxx"

namespace andres {
namespace ilp {

template<class T>
class LinearExpression {
public:
    typedef T Value;

    LinearExpression()
        :   variableIndices_(),
            coefficients_()
        {}

    template<class VariableIterator, class CoefficientIterator>
    LinearExpression(
        VariableIterator vit, 
        VariableIterator vitEnd, 
        CoefficientIterator cit
    )
        :   variableIndices_(vit, vitEnd),
            coefficients_(cit, cit + std::distance(vit, vitEnd))
        {}

    LinearExpression(
        const size_t variableIndex,
        const Value coefficient
    )
        :   variableIndices_(1, variableIndex),
            coefficients_(1, coefficient)
        {}

    template<class VariableIterator, class CoefficientIterator>
    void assign(
        VariableIterator vit, 
        VariableIterator vitEnd, 
        CoefficientIterator cit
    ) 
        { 
            variableIndices_.assign(vit, vitEnd);
            coefficients_.assign(cit, cit + std::distance(vit, vitEnd));
        }

    template<class VariableIterator, class CoefficientIterator>
    void assign(
        const size_t variableIndex,
        const Value coefficient
    )
        {
            variableIndices_.resize(1);
            variableIndices_[0] = variableIndex;
            coefficients_.resize(1);
            coefficients_[0] = coefficient;
        }

    size_t size() const
        { return variableIndices_.size(); }
    size_t index(const size_t j) const
        { return variableIndices_[j]; }
    size_t& index(const size_t j) 
        { return variableIndices_[j]; }
    Value coefficient(const size_t j) const
        { return coefficients_[j]; }
    Value& coefficient(const size_t j) 
        { return coefficients_[j]; }    

private:
    std::vector<size_t> variableIndices_;
    std::vector<Value> coefficients_;
};

template<class T>
class LinearConstraint {
public:
    typedef T Value;

    LinearConstraint()
        :   linearExpression_(),
            lowerBound_(),
            upperBound_()
        {}

    template<class VariableIterator, class CoefficientIterator>
    LinearConstraint(
        VariableIterator vit, 
        VariableIterator vitEnd, 
        CoefficientIterator cit,
        const Value lowerBound,
        const Value upperBound
    )
        :   linearExpression_(vit, vitEnd, cit),
            lowerBound_(lowerBound),
            upperBound_(upperBound)
        {}

    LinearConstraint(
        const size_t variableIndex,
        const Value lowerBound,
        const Value upperBound
    )
        :   linearExpression_(variableIndex, static_cast<Value>(1)),
            lowerBound_(lowerBound),
            upperBound_(upperBound)
        {}

    template<class VariableIterator, class CoefficientIterator>
    void assign(
        VariableIterator vit, 
        VariableIterator vitEnd, 
        CoefficientIterator cit,
        const Value lowerBound,
        const Value upperBound
    ) 
        { 
            linearExpression_.assign(vit, vitEnd, cit);
            lowerBound_ = lowerBound;
            upperBound_ = upperBound;
        }

    template<class VariableIterator, class CoefficientIterator>
    void assign(
        const size_t variableIndex,
        const Value lowerBound,
        const Value upperBound
    )
        {
            linearExpression_.assign(variableIndex, static_cast<Value>(1)),
            lowerBound_ = lowerBound;
            upperBound_ = upperBound;
        }

    size_t size() const
        { return linearExpression_.size(); }
    size_t index(const size_t j) const
        { return linearExpression_.index(j); }
    size_t& index(const size_t j) 
        { return linearExpression_.index(j); }
    Value coefficient(const size_t j) const
        { return linearExpression_.coefficient(j); }
    Value& coefficient(const size_t j) 
        { return linearExpression_.coefficient(j); }    
    Value lowerBound() const
        { return lowerBound_; }
    Value& lowerBound()
        { return lowerBound_; }
    Value upperBound() const
        { return upperBound_; }
    Value& upperBound()
        { return upperBound_; }

private:
    LinearExpression<Value> linearExpression_;
    Value lowerBound_;
    Value upperBound_;
};

template<class ILP>
class ILPCC {
public:
    typedef typename ILP::Value Value;
    enum VariableType {Int, Float, Bool};
    typedef LinearConstraint<Value> Constraint;

    struct IdleVisitor {
        void decomposition(const std::vector<std::vector<size_t> >&) const
            {}
        void beginProblem(
            const size_t, 
            const std::vector<size_t>&,
            const std::set<size_t>&
        ) const {}
        void endProblem(const size_t) const {}
    };

    struct VerboseVisitor {
        void decomposition(const std::vector<std::vector<size_t> >&) const
            {}
        void beginProblem(
            const size_t pi,
            const std::vector<size_t>& variableIndices,
            const std::set<size_t>& constraintIndices
        ) const 
            {
	      /*std::cout << "solving problem " << pi << " with " 
                    << variableIndices.size() << " variables and " 
                    << constraintIndices.size() << " constraints..." 
	            << std::endl;
	      */
            }
        void endProblem(const size_t pi) const
      { }//std::cout << "solved problem " << pi << "." << std::endl; }
    };

    // constructor
    ILPCC();
    void clear();

    // access properties
    size_t numberOfVariables() const; 
    Value value(const size_t) const;
    Value absoluteGap() const;
    Value relativeGap() const;
    size_t numberOfConstraints() const;

    // manipulate properties
    void setAbsoluteGap(const Value);
    void setRelativeGap(const Value);
    void setVerbosity(const bool);

    // set up a problem
    size_t addVariable(const VariableType); 
    size_t addVariables(const size_t, const VariableType);
    void setVariableType(const size_t, const VariableType);
    void setLinearCoefficient(const size_t, const Value);
    template<class CoefficientIterator>
        void setLinearCoefficients(CoefficientIterator);
    template<class CoefficientIterator>
        void setLinearCoefficients(CoefficientIterator, CoefficientIterator);
    template<class VariableIterator, class CoefficientIterator>
        void setLinearCoefficients(VariableIterator, VariableIterator,
            CoefficientIterator);
    void addConstraint(const size_t, const Value, const Value); 
    template<class VariableIndexIterator, class CoefficientIterator>
        void addConstraint(VariableIndexIterator, VariableIndexIterator,
                           CoefficientIterator, const Value, const Value);
    template<class VariableIndexIterator, class CoefficientIterator>
    	bool evaluateConstraint(VariableIndexIterator, VariableIndexIterator,
                           CoefficientIterator, const Value, const Value);
    // solve a problem
    template<class Iterator>
        void setStart(Iterator);    
    template<class Visitor>
        Value optimize(Visitor&);
    Value optimize();
  
private:
    struct Variable {
        Variable(const VariableType variableType = Bool)
            :   type_(variableType),
                coefficient_(),
                needsUpdate_(true),
                constraintIndices_(),
                value_(),
                valueStart_()
            {}

        VariableType type_;
        Value coefficient_;
        bool needsUpdate_;
        std::set<size_t> constraintIndices_;
        Value value_;
        Value valueStart_;
    };

    template<class VariableIndexIterator, class CoefficientIterator>
        void printConstraint(VariableIndexIterator, VariableIndexIterator,
                           CoefficientIterator, const Value, const Value) const;
    void printConstraint(const Constraint&) const;

    Value absoluteGap_;
    Value relativeGap_;
    bool verbosity_;
    std::vector<Variable> variables_;
    Partition<size_t> variablePartition_;
    std::vector<Constraint> constraints_;
};

template<class ILP>
inline 
ILPCC<ILP>::ILPCC() 
:   absoluteGap_(0),
    relativeGap_(0),
    verbosity_(false),
    variables_(),
    variablePartition_(),
    constraints_()
{}

template<class ILP>
inline void
ILPCC<ILP>::clear() {
    absoluteGap_ = 0;
    relativeGap_ = 0;
    verbosity_ = false;
    variables_.clear();
    variablePartition_.assign();
    constraints_.clear();
}

template<class ILP>
inline size_t 
ILPCC<ILP>::numberOfVariables() const {
    return variables_.size();
}

template<class ILP>
inline size_t
ILPCC<ILP>::numberOfConstraints() const {
	return constraints_.size();
}

template<class ILP>
inline typename ILPCC<ILP>::Value 
ILPCC<ILP>::value(
    const size_t j
) const {
    return variables_[j].value_;
}

template<class ILP>
inline typename ILPCC<ILP>::Value 
ILPCC<ILP>::absoluteGap() const {
    return absoluteGap_;
}

template<class ILP>
inline typename ILPCC<ILP>::Value 
ILPCC<ILP>::relativeGap() const {
    return relativeGap_;
}

template<class ILP>
inline void 
ILPCC<ILP>::setAbsoluteGap(
    const Value absoluteGap
) {
    absoluteGap_ = absoluteGap;
}

template<class ILP>
inline void 
ILPCC<ILP>::setRelativeGap(
    const Value relativeGap
) {
    relativeGap_ = relativeGap;
}

template<class ILP>
inline void 
ILPCC<ILP>::setVerbosity(
    const bool verbosity
) {
    verbosity_ = verbosity;
}

template<class ILP>
inline size_t 
ILPCC<ILP>::addVariable(
    const VariableType variableType
) {
    return addVariables(1, variableType);
}

template<class ILP>
inline size_t 
ILPCC<ILP>::addVariables(
    const size_t number, 
    const VariableType variableType
) {
    const size_t offset = numberOfVariables();
    variables_.resize(offset + number, Variable(variableType));
    variablePartition_.insert(number);
    assert(variables_.size() == variablePartition_.numberOfElements());
    return offset;
}

template<class ILP>
inline void 
ILPCC<ILP>::setVariableType(
    const size_t variableIndex, 
    const VariableType variableType
) {
    variables_[variableIndex].type_ = variableType;
    variables_[variableIndex].needsUpdate_ = true;
}

template<class ILP>
inline void 
ILPCC<ILP>::setLinearCoefficient(
    const size_t variableIndex, 
    const Value coefficient
) {
    variables_[variableIndex].coefficient_ = coefficient;
    variables_[variableIndex].needsUpdate_ = true;
}

template<class ILP>
template<class CoefficientIterator>
inline void 
ILPCC<ILP>::setLinearCoefficients(
    CoefficientIterator cit
) {
    for(size_t vi = 0; vi < numberOfVariables(); ++vi, ++cit) {
        variables_[vi].coefficient_ = *cit;
        variables_[vi].needsUpdate_ = true;
    }
}

template<class ILP>
template<class CoefficientIterator>
inline void 
ILPCC<ILP>::setLinearCoefficients(
    CoefficientIterator cit,
    CoefficientIterator citEnd
) {
    assert(std::distance(cit, citEnd) < numberOfVariables());
    for(size_t vi = 0; cit != citEnd; ++vi, ++cit) {
        variables_[vi].coefficient_ = *cit;
        variables_[vi].needsUpdate_ = true;
    }
}

template<class ILP>
template<class VariableIterator, class CoefficientIterator>
inline void 
ILPCC<ILP>::setLinearCoefficients(
    VariableIterator vit, 
    VariableIterator vitEnd,
    CoefficientIterator cit
) {
    for(; vit != vitEnd; ++vit, ++cit) {
        variables_[*vit].coefficient_ = *cit;
        variables_[*vit].needsUpdate_ = true;
    }
}

template<class ILP>
inline void 
ILPCC<ILP>::addConstraint(
    const size_t variableIndex,  
    const Value lowerBound, 
    const Value upperBound
) {
    constraints_.push_back(
        Constraint(variableIndex, lowerBound, upperBound)
    );
    const size_t constraintIndex = constraints_.size() - 1;
    variables_[variableIndex].constraintIndices_.insert(constraintIndex);
    variables_[variableIndex].needsUpdate_ = true;

    // debug output
    /*
    std::cout << "adding constraint: " << lowerBound 
        << " <= x" << variableIndex
        << " <= " << upperBound << std::endl;
//    */
}

template<class ILP>
template<class VariableIndexIterator, class CoefficientIterator>
inline void 
ILPCC<ILP>::addConstraint(
    VariableIndexIterator vit, 
    VariableIndexIterator vitEnd,
    CoefficientIterator cit, 
    const Value lowerBound, 
    const Value upperBound
) {
    // begin debug output
	/*
		std::cout << "adding constraint: ";
		printConstraint(vit, vitEnd, cit, lowerBound, upperBound);
//	*/
    // end debug output

    constraints_.push_back(
        Constraint(vit, vitEnd, cit, lowerBound, upperBound)
    );
    const size_t constraintIndex = constraints_.size() - 1;
    for(; vit != vitEnd; ++vit) {
        variables_[*vit].constraintIndices_.insert(constraintIndex);
        variables_[*vit].needsUpdate_ = true;
        const VariableIndexIterator nextVit = vit + 1;
        if(nextVit != vitEnd) {
            variablePartition_.merge(*vit, *nextVit);
        }
    }
}

template<class ILP>
template<class Iterator>
inline void 
ILPCC<ILP>::setStart(
    Iterator valueIterator
) {
    for(size_t j = 0; j < numberOfVariables(); ++j, ++valueIterator) {
        variables_[j].valueStart_ = *valueIterator;
    }
}

template<class ILP>
inline typename ILPCC<ILP>::Value 
ILPCC<ILP>::optimize() {
    VerboseVisitor visitor;
    return optimize<VerboseVisitor>(visitor);
}

template<class ILP>
template<class Visitor>
inline typename ILPCC<ILP>::Value 
ILPCC<ILP>::optimize(
    Visitor& visitor
) {
    std::vector<size_t> problemsOfVariables(variablePartition_.numberOfElements());
    variablePartition_.elementLabeling(problemsOfVariables.begin());

    // - make variable partition explicit
    // - determine which sub-problems need to be re-solved
    const size_t numberOfProblems = variablePartition_.numberOfSets();
    std::vector<std::vector<size_t> > variablesOfProblems(numberOfProblems);
    std::vector<bool> updateRequirementOfProblems(numberOfProblems);
    for(size_t vi = 0; vi < numberOfVariables(); ++vi) {
        const size_t pi = problemsOfVariables[vi];
        variablesOfProblems[pi].push_back(vi);
        if(variables_[vi].needsUpdate_) {
            updateRequirementOfProblems[pi] = true;
        }
    }
    // debug output
//    std::cout << ">> #Problems: " << numberOfProblems << std::endl;
//    std::cout << "(max. size: " << << ")" << std::endl;

    visitor.decomposition(variablesOfProblems);
    Value objectiveValue = 0;

    #pragma omp parallel 
    {
        std::vector<size_t> viSubOfVi(numberOfVariables()); 
        // std::map<size_t, size_t> viSubOfVi; 
        std::set<size_t> constraintIndicesSub;  
        std::vector<Value> coefficientsSub;  
        coefficientsSub.reserve(10);
        #pragma omp for schedule(dynamic)
        // for(ptrdiff_t pi = 0; pi < static_cast<ptrdiff_t>(numberOfProblems); ++pi) { // for every sub-problem
        for(std::ptrdiff_t pi = 0; pi < static_cast<std::ptrdiff_t>(numberOfProblems); ++pi) { // for every sub-problem
//             std::cout << "sub-problem " << pi << std::endl; // debug output

            if(!updateRequirementOfProblems[pi]) { 
//                 std::cout << "already up to date." << std::endl; // debug output
                continue;
            }
//             std::cout << "needs to be solved." << std::endl; // debug output
  
            // - make explicit the mapping from variable indices to sub-problem variable indices
            // - make explicit the indices of constraints on the sub-problem
            const std::vector<size_t>& variablesSub = variablesOfProblems[pi];
            const VariableType variableTypeSub = variables_[variablesSub[0]].type_;
            bool uniqueVariableType = true;
        
            // viSubOfVi.clear(); // for a map
            constraintIndicesSub.clear();
            for(size_t viSub = 0; viSub < variablesSub.size(); ++viSub) {
                const size_t vi = variablesSub[viSub];
                viSubOfVi[vi] = viSub;
                if(variables_[vi].type_ != variableTypeSub) {
                    uniqueVariableType = false;
                }
                for(std::set<size_t>::const_iterator cit = variables_[vi].constraintIndices_.begin(); 
                cit != variables_[vi].constraintIndices_.end(); ++cit) {
                    constraintIndicesSub.insert(*cit);
                }
            }

            // begin debug output
            /*
            std::cout << "variable indices:";
            for(size_t j = 0; j < variablesSub.size(); ++j) {
                std::cout << ' ' << variablesSub[j] << '(' << j << ')';
            }
            std::cout << std::endl << "constraints:" << std::endl;
            for(std::set<size_t>::const_iterator cit = constraintIndicesSub.begin();
            cit != constraintIndicesSub.end(); ++cit) {
                const Constraint& constraint = constraints_[*cit];
                std::cout << *cit << ": ";
                printConstraint(constraints_[*cit]);
            }
            */
            // end debug output

            #pragma omp critical
            visitor.beginProblem(pi, variablesSub, constraintIndicesSub);

            if(variablesSub.size() == 1 && variables_[variablesSub[0]].type_ == Bool) { 
                // determine feasibility of solutions 0 and 1
                bool feasible[] = {true, true};
                for(std::set<size_t>::const_iterator cit = constraintIndicesSub.begin();
                cit != constraintIndicesSub.end(); ++cit) {
                    const Constraint& constraint = constraints_[*cit];
                    assert(constraint.size() == 1);

                    // evaluate lower bound
                    if(constraint.lowerBound() > constraint.coefficient(0) * 0) {
                        feasible[0] = false;
                    }               
                    if(constraint.lowerBound() > constraint.coefficient(0) * 1) {
                        feasible[1] = false;
                    }

                    // evaluate upper bound
                    if(constraint.coefficient(0) * 0 > constraint.upperBound()) {
                        feasible[0] = false;
                    }                
                    if(constraint.coefficient(0) * 1 > constraint.upperBound()) {
                        feasible[1] = false;
                    }
                }

                // determine optimal feasible solution
                Value solutionSub;
                Value objectiveValueSub;
                const Value coefficient = variables_[variablesSub[0]].coefficient_;
                if(feasible[1]) {
                    if(feasible[0]) {
                        if(0 <= coefficient) { // if 0 is better than 1
                            solutionSub = 0;
                            objectiveValueSub = 0;
                        }
                        else {
                            solutionSub = 1;
                            objectiveValueSub = coefficient;
                        }
                    }
                    else {
                        solutionSub = 1;
                        objectiveValueSub = coefficient;
                    }
                }
                else {
                    if(feasible[0]) {
                        solutionSub = 0;
                        objectiveValueSub = 0;
                    }
                    else {
                        throw std::runtime_error("infeasible.");
                    }
                }

                #pragma omp atomic
                objectiveValue += objectiveValueSub;

                // store solution of sub-problem and flag variable as up-to-date
                const size_t vi = variablesSub[0];
                variables_[vi].value_ = solutionSub;
                variables_[vi].needsUpdate_ = false;
            }
            else if(variablesSub.size() == 1 && variables_[variablesSub[0]].type_ == Float) { 
                // determine bounds on the variable
                Value smallestUpperBound = std::numeric_limits<Value>::infinity();
                Value greatestLowerBound = -std::numeric_limits<Value>::infinity();
                for(std::set<size_t>::const_iterator cit = constraintIndicesSub.begin();
                cit != constraintIndicesSub.end(); ++cit) {
                    const Constraint& constraint = constraints_[*cit];
                    assert(constraint.size() == 1);
                    if(constraint.coefficient(0) > 0) {
                        const double upperBound = constraint.upperBound() / constraint.coefficient(0);
                        if(upperBound < smallestUpperBound) {
                            smallestUpperBound = upperBound;
                        }
                        const double lowerBound = constraint.lowerBound() / constraint.coefficient(0);
                        if(lowerBound > greatestLowerBound) {
                            greatestLowerBound = lowerBound;
                        }
                    }
                    else if(constraint.coefficient(0) < 0) {
                        const double upperBound = constraint.lowerBound() / constraint.coefficient(0);
                        if(upperBound < smallestUpperBound) {
                            smallestUpperBound = upperBound;
                        }
                        const double lowerBound = constraint.upperBound() / constraint.coefficient(0);
                        if(lowerBound > greatestLowerBound) {
                            greatestLowerBound = lowerBound;
                        }
                    }
                    else { // constraint.coefficient() == 0
                        if(constraint.lowerBound() > 0 && constraint.upperBound() < 0) {
                            throw std::runtime_error("infeasbile.");
                        }
                    }
                }

                // solve problem
                if(smallestUpperBound == std::numeric_limits<Value>::infinity()
                || greatestLowerBound == -std::numeric_limits<Value>::infinity()) {
                    throw std::runtime_error("unbounded.");
                }
                else if(greatestLowerBound > smallestUpperBound) {
                    throw std::runtime_error("infeasible.");
                }
                else {
                    // store solution of sub-problem and flag variable as up-to-date
                    assert(greatestLowerBound < smallestUpperBound);
                    const size_t vi = variablesSub[0];
                    variables_[vi].needsUpdate_ = false;
                    const Value coefficient = variables_[variablesSub[0]].coefficient_;
                    if(coefficient > 0) {
                        variables_[vi].value_ = greatestLowerBound;
                    }
                    else if(coefficient < 0) {
                        variables_[vi].value_ = smallestUpperBound;
                    }
                    else { // coefficient == 0
                        variables_[vi].value_ = (smallestUpperBound - greatestLowerBound) / 2;
                    }
                }
            }
            else { 
                // initialize sub-problem and solver
                ILP ilp;
                ilp.setNumberOfThreads(4);
                ilp.setAbsoluteGap(absoluteGap_);
                ilp.setRelativeGap(relativeGap_);

                // add sub-problem variables to solver 
                if(uniqueVariableType) {
                    switch(variableTypeSub) {
                    case Float:
                        ilp.addVariables(variablesSub.size(), ILP::Float);
                        break;
                    case Int:
                        ilp.addVariables(variablesSub.size(), ILP::Int);
                        break;
                    case Bool:
                        ilp.addVariables(variablesSub.size(), ILP::Bool);
                        break;
                    default:
                        break;
                    }
                }
                else {
                    for(size_t viSub = 0; viSub < variablesSub.size(); ++viSub) {
                        const size_t vi = variablesSub[viSub];
                        switch(variables_[vi].type_) {
                        case Float:
                            ilp.addVariable(ILP::Float);
                            break;
                        case Int:
                            ilp.addVariable(ILP::Int);
                            break;
                        case Bool:
                            ilp.addVariable(ILP::Bool);
                            break;
                        default:
                            break;
                        }
                    }
                }

                // set objective function and warm-start for sub-problem
                {
                    coefficientsSub.resize(variablesSub.size());
                    for(size_t viSub = 0; viSub < variablesSub.size(); ++viSub) {
                        const size_t vi = variablesSub[viSub];
                        coefficientsSub[viSub] = variables_[vi].coefficient_;
                        ilp.registerStart(viSub, variables_[vi].valueStart_);
                        ilp.registerStart(viSub, 0.0);
                        // debug output
//                        std::cout << vi << ":" << variables_[vi].valueStart_ << "/" << variables_[vi].value_ <<"/" << variables_[vi].coefficient_ << ", " ;	// debug
                    }
//                    std::cout << std::endl;		// debug
                    ilp.setLinearCoefficients(coefficientsSub.begin());
                }

                // add constraints for sub-problem
                std::vector<size_t> variableIndices;    
                variableIndices.reserve(30);
                std::vector<Value> coefficients;        
                coefficients.reserve(30);
                for(std::set<size_t>::const_iterator cit = constraintIndicesSub.begin();
                cit != constraintIndicesSub.end(); ++cit) {
                    const Constraint& constraint = constraints_[*cit];
                    variableIndices.resize(constraint.size());
                    coefficients.resize(constraint.size());
                    for(size_t j = 0; j < constraint.size(); ++j) {
                        const size_t vi = constraint.index(j);
                        assert(problemsOfVariables[vi] == pi);
                        variableIndices[j] = viSubOfVi[vi];
                        coefficients[j] = constraint.coefficient(j);
                    }
                    ilp.addConstraint(
                        variableIndices.begin(), variableIndices.end(),
                        coefficients.begin(), 
                        constraint.lowerBound(), constraint.upperBound()
                    );
                }

                // solve sub-problem
                const size_t objectiveValueSub = ilp.optimize();

                #pragma omp atomic
                objectiveValue += objectiveValueSub;

                // store solution of sub-problem and flag variables as up-to-date
                for(size_t viSub = 0; viSub < variablesSub.size(); ++viSub) {
                    const size_t vi = variablesSub[viSub];
                    variables_[vi].value_ = ilp.value(viSub);
                    variables_[vi].needsUpdate_ = false;
                }
            }

            #pragma omp critical
            visitor.endProblem(pi);
        }
    }
    // FIXME
    return objectiveValue;
}

template<class ILP>
template<class VariableIndexIterator, class CoefficientIterator>
void 
ILPCC<ILP>::printConstraint(
    VariableIndexIterator vit, 
    VariableIndexIterator vitEnd,
    CoefficientIterator cit, 
    const Value lowerBound, 
    const Value upperBound
) const {
    std::cout << lowerBound << " <=";
    for(; vit != vitEnd; ++vit, ++cit) {
        if(*cit < 0) {
            if(*cit == -1) {
                std::cout << " -";
            }
            else {
                std::cout << " - " << -(*cit);
            }
        }
        else if(*cit > 0) {
            if(*cit == 1) {
                std::cout << " +";
            }
            else {
                std::cout << " + " << *cit;
            }
        }
        std::cout << " x" << *vit;
    }
    std::cout << " <= " << upperBound << std::endl;
}

template<class ILP>
void 
ILPCC<ILP>::printConstraint(
    const Constraint& constraint
) const {
    std::cout << constraint.lowerBound() << " <=";
    for(size_t j = 0; j < constraint.size(); ++j) {
        if(constraint.coefficient(j) < 0) {
            if(constraint.coefficient(j) == -1) {
                std::cout << " -";
            }
            else {
                std::cout << " - " << -constraint.coefficient(j);
            }
        }
        else if(constraint.coefficient(j) > 0) {
            if(constraint.coefficient(j) == 1) {
                std::cout << " +";
            }
            else {
                std::cout << " + " << constraint.coefficient(j);
            }
        }
        std::cout << " x" << constraint.index(j);
    }
    std::cout << " <= " << constraint.upperBound() << std::endl;
}

} // namespace ilp
} // namespace andres

#endif // #ifndef ANDRES_ILP_ILPCC_HXX
