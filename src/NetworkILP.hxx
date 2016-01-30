#ifndef NETWORK_ILP_HXX
#define NETWORK_ILP_HXX

// std
#include <iostream>
#include <deque>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>

// graph
#include "andres/graph/digraph.hxx"
#include "andres/graph/components.hxx"
#include "andres/graph/shortest-paths.hxx"

namespace netopt {
namespace ilp {

template<class ILP>
class NetworkILP {
public:

	/// helper struct for connected-components analysis.
	template<class T>
	struct Mask {
		const T * m_weights;
		explicit Mask(const T * weights) {
			m_weights = weights;
		}
		bool vertex(const size_t v) const {
			return m_weights->at(v) < 0;
		}
		bool edge(const size_t e) const {
			return true;
		}
	};

	struct AuxTuple {
		size_t m_aux;
		std::vector<size_t> m_members;
		AuxTuple(
			size_t aux,
			std::vector<size_t> members
		) {
			m_aux = aux;
			m_members = members;
		}
		// utility for debugging.
		void describe() const {
			std::cout << "aux=" << m_aux << " members={";
			for (auto member : m_members) {
				std::cout << member << ", ";
			}
			std::cout << "}" << std::endl;
		}
	};

	// typedefs
	typedef std::vector<size_t> NodeVectorType;
	typedef std::vector<typename ILP::Value> ValueVectorType;
	typedef andres::graph::Digraph<> GraphType;
	typedef andres::graph::ComponentsByPartition<GraphType> ComponentsType;
	typedef std::deque<size_t> IndicatorVectorType;
	typedef std::deque<typename ILP::Value> WeightType;
	typedef std::set<size_t> NodeSetType;
	typedef std::deque<AuxTuple> AuxTupleDequeType;

	// constructor/destructor.
	NetworkILP() {}
	explicit NetworkILP(size_t n);
	~NetworkILP() {} 	// empty.

	// problem construction.
	void addNode(float);
	void addEdge(size_t, size_t);
	void addNodes(const ValueVectorType &);
	void addEdges(const NodeVectorType &, const NodeVectorType &);
	void addDirectedEdges(const NodeVectorType &, const NodeVectorType &);
	void addPairwiseWeight(size_t, size_t, typename ILP::Value);
	void addPairwiseWeights(const NodeVectorType &, const NodeVectorType &, const ValueVectorType &);
	void addTripletWeight(size_t, size_t, size_t, typename ILP::Value);
	void addTripletWeights(const NodeVectorType &, const NodeVectorType &,
			const NodeVectorType &, const ValueVectorType &);
	void update();

	// constraints.
	void activateLazyConnectivityConstraints();	// this only sets a flag.
	void addCircleConstraint(const NodeVectorType &);

	// find a solution.
	typename ILP::Value optimize();
	typename ILP::Value optimize(const bool);

	// access solution.
	ValueVectorType getValues() const;
	ValueVectorType getAllValues() const;

	// simple getters/setters.
	void setAbsoluteGap(typename ILP::Value);
	void setRelativeGap(const typename ILP::Value);
	typename ILP::Value getAbsoluteGap() const;
	typename ILP::Value getRelativeGap() const;
	size_t getNumberOfNodes() const;
	size_t getNumberOfEdges() const;
	size_t getNumberOfVariables() const;
	size_t getNumberOfConstraints() const;

	void setRoot(size_t);
	size_t getRoot();

private:
	typename ILP::Value lazyNodeConnectivityConstraints();
	void addSetLazyNodeConstraint(const NodeSetType & , const NodeSetType &);
	void addSingleLazyNodeConstraint(size_t, const NodeSetType &);
	void determineBestRoot();
	void addSmallerThanConstraint(const size_t, const size_t,
			const typename ILP::Value, const typename ILP::Value);
	void addPairAuxiliaryConstraints(const size_t, const size_t, const size_t);
	void addAuxiliaryNANDConstraints(const size_t, const size_t, const size_t);
	void addMutuallyExclusiveConstraint(const size_t, const size_t);
	bool checkAuxiliaryConstraints();
	void addAuxiliaryConstraintsFromTuple(const AuxTuple &);

	WeightType m_weights;
	ILP m_ilp;
	GraphType m_graph;
	bool m_isUpdated{false};
	IndicatorVectorType m_indicators;
	AuxTupleDequeType m_unconstrained;
	size_t m_root{0};
	bool m_relaxation{false};
	bool m_lazy{false};
};


template<class ILP>
inline
NetworkILP<ILP>::NetworkILP(size_t n) {
	m_weights = WeightType(n, 0);
}

/// adds a decision node to the problem.
///
/// Note: A decision node may, for example correspond to a vessel segment.
template<class ILP>
inline void
NetworkILP<ILP>::addNode(float weight) {
	size_t id = m_graph.insertVertex();
	if (id < m_weights.size()) {
		m_weights[id] = weight;
	} else {
		m_weights.push_back(weight);
	}
	m_isUpdated = false;
}

/// adds an edge between two decision nodes.
///
/// These are used to determine neighbouring segments.
template<class ILP>
inline void NetworkILP<ILP>::addEdge(size_t id1, size_t id2) {
	m_graph.insertEdge(id1, id2);
	m_isUpdated = false;
}

/// adds a set of nodes with their weights at once.
///
template<class ILP>
inline void
NetworkILP<ILP>::addNodes(
		const ValueVectorType & weights
) {
	for (auto w : weights) {
		this->addNode(w);
	}
}

/// add set of edges.
/// note that all edge pairs are added bidirectionally.
///
template<class ILP>
inline void
NetworkILP<ILP>::addEdges(
		const NodeVectorType & x,
		const NodeVectorType & y
) {
	if (x.size() != y.size()) {
		throw std::invalid_argument(
				"start and stop node ids must contain the same number of elements");
	}
	m_graph.reserveEdges(x.size() * 2);
	for (size_t ii = 0; ii < x.size(); ii++) {
		this->addEdge(x[ii], y[ii]);
		this->addEdge(y[ii], x[ii]);
	}
}

/// add directed edges.
///
template<class ILP>
inline void
NetworkILP<ILP>::addDirectedEdges(
		const NodeVectorType & x,
		const NodeVectorType & y
) {
	if (x.size() != y.size()) {
		throw std::invalid_argument(
				"start and stop node ids must contain the same number of elements");
	}
	m_graph.reserveEdges(x.size());
	for (size_t ii = 0; ii < x.size(); ii++) {
			this->addEdge(x[ii], y[ii]);
	}
}

/// adds a pairwise auxiliary variable with given weight for
/// two given node variables. The auxiliary variable z is constrained
/// to
/// 	z = x_i * x_j
///
/// NOTE: update has to be called before doing this!
template<class ILP>
inline void
NetworkILP<ILP>::addPairwiseWeight(
		size_t x,
		size_t y,
		typename ILP::Value weight
) {
	// during the optimization, we check if the constraints are necessary.
	// and will call addPairAuxiliaryConstraints()
	// if necessary.
	size_t aux = m_ilp.addVariable(ILP::Bool);
	m_ilp.setLinearCoefficient(aux, weight);
	m_unconstrained.push_back(AuxTuple(aux, {m_indicators[x], m_indicators[y]}));
}


/// adds a triplet auxiliary variable with given weight for
/// three node variables. The auxiliary variable z is constrained
/// to
/// 	z = x_i * x_j * x_k
///
/// NOTE: update has to be called before doing this!
template<class ILP>
inline void
NetworkILP<ILP>::addTripletWeight(
		size_t x,
		size_t y,
		size_t z,
		typename ILP::Value weight
) {
	size_t aux = m_ilp.addVariable(ILP::Bool);
	m_ilp.setLinearCoefficient(aux, weight);
	m_unconstrained.push_back(AuxTuple(aux,
			{m_indicators[x], m_indicators[y], m_indicators[z]}));
}

/// vectorized interface to add pairwise weights.
///
template<class ILP>
inline void
NetworkILP<ILP>::addPairwiseWeights(
		const NodeVectorType & xi,
		const NodeVectorType & xj,
		const ValueVectorType & w
) {
	if (xi.size() != xj.size() || w.size() != xi.size()) {
		throw std::invalid_argument(
				"all input vectors must contain the same number of elements");
	}
	for (size_t cc = 0; cc < xi.size(); cc++) {
		this->addPairwiseWeight(xi[cc], xj[cc], w[cc]);
	}
}

/// vectorized interface to add weights for triplets.
///
template<class ILP>
inline void
NetworkILP<ILP>::addTripletWeights(
		const NodeVectorType & xi,
		const NodeVectorType & xj,
		const NodeVectorType & xk,
		const ValueVectorType & w
) {
	if (xi.size() != xj.size() || xi.size() != xk.size()
			|| xi.size() != w.size()) {
			throw std::invalid_argument(
				"all input vectors must contain the same number of elements");
	}
	for (size_t cc = 0; cc < xi.size(); cc++) {
		this->addTripletWeight(xi[cc], xj[cc], xk[cc], w[cc]);
	}
}

/// adds missing variables to the ILP and updates _all_ weights.
///
template<class ILP>
void
NetworkILP<ILP>::update() {
	if (m_isUpdated) {
		return;
	}
	const size_t missing = m_weights.size() - m_indicators.size();
	if (missing > 0) {
		size_t baseId;
		if (m_relaxation) {
			baseId = m_ilp.addVariables(missing, ILP::Float);
			// limit relaxed indicators to [0,1]
			for (int ii = 0; ii < missing; ii++) {
				m_ilp.addConstraint(ii + baseId, 0, 1);
			}
		} else {
			baseId = m_ilp.addVariables(missing, ILP::Bool);
		}
		// we keep track of the ids of indicator variables here.
		for (int ii = 0; ii < missing; ii++) {
			m_indicators.push_back(ii+baseId);
		}
	}

	// if we have no auxiliary variables, then we can call the efficient
	// setter function.
	auto wit = m_weights.begin();
	for (auto it = m_indicators.begin();
			it != m_indicators.end() && wit != m_weights.end(); ++it, ++wit) {
		m_ilp.setLinearCoefficient(*it, *wit);
	}
	m_isUpdated = true;
}


/// Selects the node with the lowest weight in the largest connected
/// component as root node.
///
template<class ILP>
inline void
NetworkILP<ILP>::determineBestRoot() {
	// Connected components analysis.
	ComponentsType components;
	components.build(m_graph, Mask<WeightType>(&m_weights));
	auto labelling = NodeVectorType(m_graph.numberOfVertices());
	components.partition_.elementLabeling(labelling.begin());

	size_t maxMembers = 0;
	for (size_t nn = 0; nn < components.partition_.numberOfSets(); nn++) {
		const size_t n_members = std::count(labelling.begin(), labelling.end(), nn);
		if (maxMembers > n_members) {
			continue;
		} else {
			maxMembers = n_members;
			size_t ii = 0;
			typename ILP::Value minWeight = 0;
			for (auto it = m_weights.begin(); it != m_weights.end(); ++it, ++ii) {
				if ((labelling[ii] == nn) && (*it < minWeight)) {
					m_root = ii;
					minWeight = *it;
				}
			}
		}
	}
}


/// adds x_ii <= x_jj constraint
///
template<class ILP>
inline void
NetworkILP<ILP>::addSmallerThanConstraint(
		const size_t ii,
		const size_t jj,
		const typename ILP::Value lower,
		const typename ILP::Value upper
) {
	typename ILP::Value coeffs[] = {1, -1};
	size_t variableIndices[] = {jj, ii};
	m_ilp.addConstraint(variableIndices, variableIndices+2,
			coeffs, lower, upper);
}

/// adds a mutually-exclusive constraint
///
template<class ILP>
inline void
NetworkILP<ILP>::addMutuallyExclusiveConstraint(
		const size_t ii,
		const size_t jj
) {
	typename ILP::Value coeffs[] = {1, 1};
	typename ILP::Value lower = 0;
	typename ILP::Value upper = 1;
	size_t variableIndices[] = {ii, jj};
	m_ilp.addConstraint(variableIndices, variableIndices+2,
			coeffs, lower, upper);
}

/// manual addition of pairwise auxiliary constraints s.t.
///   y_aux == x_ii * x_jj.
///
template<class ILP>
inline void
NetworkILP<ILP>::addPairAuxiliaryConstraints(
		const size_t aux,
		const size_t ii,
		const size_t jj
) {
	// y_ij <= x_i and y_ij <= x_j
	const typename ILP::Value lower = 0;
	const typename ILP::Value upper = 1;
	this->addSmallerThanConstraint(aux, ii, lower, upper);
	this->addSmallerThanConstraint(aux, jj, lower, upper);

	// 0 <= x_i + x_j - y_ij <= 1
	typename ILP::Value coeffs[] = {-1, 1, 1};
	size_t variableIndices[] = {aux, ii, jj};
	m_ilp.addConstraint(variableIndices, variableIndices+3,
			coeffs, lower, upper);
}

/// adds auxiliary constraints for variables represented in auxtuple.
/// 	y_aux == \prod \sum x_member
///
template<class ILP>
inline void
NetworkILP<ILP>::addAuxiliaryConstraintsFromTuple(
		const AuxTuple & auxtpl
) {
	const typename ILP::Value lower = 0;
	const typename ILP::Value upper = 1;
	size_t len = auxtpl.m_members.size();
	typename ILP::Value coeffs[1 + len];
	size_t variableIndices[1 + len];
	coeffs[0] = -1;
	variableIndices[0] = auxtpl.m_aux;

	// y <= x \forall x \in {i,j, ...}
	int ii = 1;
	for (auto mit : auxtpl.m_members) {
		this->addSmallerThanConstraint(auxtpl.m_aux, mit, lower, upper);
		variableIndices[ii] = mit;
		coeffs[ii] = 1;
		ii++;
	}
	// 0 <= -y + \sum x  <= N - 1
	m_ilp.addConstraint(variableIndices, variableIndices + len + 1,
			coeffs, lower, len - 1);
}

/// searches and adds violated auxiliary constraints.
///
template<class ILP>
inline bool
NetworkILP<ILP>::checkAuxiliaryConstraints() {
	bool addedConstraint = false;
	typename ILP::Value diff;
	typename ILP::Value prod;
	// check y == \prod x for all entries.
	for (auto auxit = m_unconstrained.begin(); auxit != m_unconstrained.end(); ) {
		prod = 1;
		for (auto m : auxit->m_members) {
			prod *= m_ilp.value(m);
		}
		diff = m_ilp.value(auxit->m_aux) - prod;
		if (diff * diff >= 0.01) {	// >= eps check due to integrality slack.
			addAuxiliaryConstraintsFromTuple(*auxit);
			auxit = m_unconstrained.erase(auxit);
			addedConstraint = true;
		} else {
			++auxit;
		}
	}
	return addedConstraint;
}


/// adds the necessary linear inequality constraints for y_ij = (1-x_i)(1-x_j)
///
template<class ILP>
inline void
NetworkILP<ILP>::addAuxiliaryNANDConstraints(
		const size_t aux,
		const size_t ii,
		const size_t jj
) {
	// y_ij <= 1-x_i and y_ij <= 1-x_j
	this->addMutuallyExclusiveConstraint(aux, ii);
	this->addMutuallyExclusiveConstraint(aux, jj);

	// 1 <= x_i + x_j + y_ij <= 2
	const typename ILP::Value lower = 1;
	const typename ILP::Value upper = 2;
	typename ILP::Value coeffs[] = {1, 1, 1};
	size_t variableIndices[] = {aux, ii, jj};
	m_ilp.addConstraint(variableIndices, variableIndices+3,
			coeffs, lower, upper);
}

/// activates generation of connectivity constraints in optimize().
///
template<class ILP>
inline void
NetworkILP<ILP>::activateLazyConnectivityConstraints() {
	m_lazy = true;
}


/// Enforces connectivity by iteratively adding violated constraints.
///
template<class ILP>
inline typename ILP::Value
NetworkILP<ILP>::lazyNodeConnectivityConstraints() {
	this->update();

	// determine root node.
	if (m_root == 0) {
		this->determineBestRoot();
	}

	// prepare mask and components.
	WeightType mask(m_indicators.size(), 1);
	ComponentsType components;
	size_t nComponents;
	size_t activeComponents = 0;
	typename ILP::Value energy;

	do {
		// 1. search violated auxiliary constraints.
		do {
			energy = m_ilp.optimize();
		} while (checkAuxiliaryConstraints());

		// 2. search violated connectivity constraints.
		// update mask.
		for (auto ii = 0; ii < m_indicators.size(); ii++) {
			// we test against 0.9 in order to avoid false branching
			// due to integrality slack in the ILP solver.
			if (m_ilp.value(m_indicators[ii]) >= 0.9) {
				mask[ii] = -1;
			} else {
				mask[ii] = 1;
			}
		}

		// determine violated connectivity constraints by analysing
		// the connected components under the masked graph.
		nComponents = components.build(m_graph, Mask<WeightType>(&mask));
		NodeVectorType rep(static_cast<size_t>(nComponents));
		components.partition_.representatives(rep.begin());

		activeComponents = 0;
		for (auto repit = rep.begin(); repit != rep.end(); ++repit) {
			// if it's an inactive component or already containing the root, we
			// can skip this set.
			if ((mask[*repit] >= 0) || components.areConnected(m_root, *repit)) {
				continue;
			} else {
				activeComponents++;
			}

			// explore the neighbourhood of each representative to obtain the
			// set of both members and neighbours.
			IndicatorVectorType membersQueue;
			NodeSetType members;
			NodeSetType neighbours;
			membersQueue.push_back(*repit);
			members.insert(*repit);

			while (!membersQueue.empty()) {
				const size_t member = membersQueue.front();
				for (auto nit = m_graph.verticesToVertexBegin(member);
						nit != m_graph.verticesToVertexEnd(member); ++nit) {
					// if the adjacent node has the same partition label as the
					// current member, then it is a member itself. otherwise,
					// it has to be a neighbour.
					if (components.areConnected(member, *nit)) {
						const bool res = members.insert(*nit).second;
						if (res) {
							membersQueue.push_back(*nit);
						}
					} else {
						neighbours.insert(*nit);
					}
				}
				membersQueue.pop_front();
			}
			// add constraints for all members.
			for (auto member : members) {
				addSingleLazyNodeConstraint(member, neighbours);
			}
		}
	} while (activeComponents > 0);
	return energy;
}


/// add constraint for active neighbourhood on a single node.
///
///		0 <= \sum_j \in {\Union N(k)\{x_i}}  - x_i <= |N|
///
template<class ILP>
inline void
NetworkILP<ILP>::addSingleLazyNodeConstraint(
		size_t member,
		const NodeSetType & neighbours
) {
	typename ILP::Value lower = 0;
	typename ILP::Value upper = neighbours.size();
	const size_t n = 1 + neighbours.size();

	size_t variableIndices[n];
	typename ILP::Value coeffs[n];
	unsigned int counter = 0;
	for (auto nit = neighbours.begin();
			nit != neighbours.end(); ++nit, ++counter) {
		coeffs[counter] = 1;
		variableIndices[counter] = m_indicators[*nit];
	}
	coeffs[counter] = -1;
	variableIndices[counter] = m_indicators[member];
	m_ilp.addConstraint(variableIndices, variableIndices + n,
			coeffs, lower, upper);
}

/// prevents circle on the given variables.
///
///  \sum_{x_i \in C} x_i <= |C|-1
///
template<class ILP>
inline void
NetworkILP<ILP>::addCircleConstraint(const NodeVectorType & ids) {
	const size_t n = ids.size();
	typename ILP::Value lower = 0;
	typename ILP::Value upper = ids.size() - 1;

	size_t variableIndices[n];
	typename ILP::Value coeffs[n];
	unsigned int counter = 0;
	for (auto id : ids) {
		variableIndices[counter] = m_indicators[id];
		coeffs[counter] = 1;
		counter++;
	}
	m_ilp.addConstraint(variableIndices, variableIndices + n,
			coeffs, lower, upper);
}

/// optimize the model.
///
template<class ILP>
inline typename ILP::Value
NetworkILP<ILP>::optimize() {
	if (m_lazy) {
		return lazyNodeConnectivityConstraints();
	} else {
		typename ILP::Value val;
		do {
			val = m_ilp.optimize();
		} while (checkAuxiliaryConstraints());
		return val;
	}
}

/// optimize with verbose output.
///
template<class ILP>
inline typename ILP::Value
NetworkILP<ILP>::optimize(const bool verbose) {
	m_ilp.setVerbosity(verbose);
	return this->optimize();
}

/// get current solution.
///
template<class ILP>
inline
typename NetworkILP<ILP>::ValueVectorType
NetworkILP<ILP>::getValues() const {
	NetworkILP<ILP>::ValueVectorType values;
	values.reserve(m_indicators.size());
	for (auto ii : m_indicators) {
		values.push_back(m_ilp.value(ii));
	}
	return values;
}

/// get current solution and values of auxiliaries.
///
template<class ILP>
inline
typename NetworkILP<ILP>::ValueVectorType
NetworkILP<ILP>::getAllValues() const {
	NetworkILP<ILP>::ValueVectorType values;
	values.reserve(m_ilp.numberOfVariables());
	for (size_t ii = 0; ii < m_ilp.numberOfVariables(); ii++) {
		values.push_back(m_ilp.value(ii));
	}
	return values;
}

/// Sets the absolute gap used as termination criterion.
///
template<class ILP>
void
NetworkILP<ILP>::setAbsoluteGap(
		typename ILP::Value gap
) {
	m_ilp.setAbsoluteGap(gap);
}

/// Sets the relative gap used as termination criterion.
///
template<class ILP>
void
NetworkILP<ILP>::setRelativeGap(
		const typename ILP::Value gap
		) {
	m_ilp.setRelativeGap(gap);
}

template<class ILP>
typename ILP::Value
NetworkILP<ILP>::getAbsoluteGap() const {
	return m_ilp.absoluteGap();
}

template<class ILP>
typename ILP::Value
NetworkILP<ILP>::getRelativeGap() const {
	return m_ilp.relativeGap();
}

/// returns the number of nodes in the graph.
///
template<class ILP>
inline size_t
NetworkILP<ILP>::getNumberOfNodes() const {
	return m_graph.numberOfVertices();
}

/// returns the number of edges in the graph.
///
template<class ILP>
inline size_t
NetworkILP<ILP>::getNumberOfEdges() const {
	return m_graph.numberOfEdges();
}

/// returns the number of variables in the ilp.
///
template<class ILP>
inline size_t
NetworkILP<ILP>::getNumberOfVariables() const {
	return m_ilp.numberOfVariables();
}

/// returns the number of constraints.
///
template<class ILP>
inline size_t
NetworkILP<ILP>::getNumberOfConstraints() const {
	return m_ilp.numberOfConstraints();
}

/// sets node with given id as root node (used for connectivity).
///
template<class ILP>
inline void
NetworkILP<ILP>::setRoot(size_t id) {
	if (m_graph.numberOfVertices() <= id) {
		throw std::invalid_argument(
				"root node id must be existing!");
	}
	m_root = id;
}

/// return the current root node.
///
template<class ILP>
inline size_t
NetworkILP<ILP>::getRoot() {
	return m_root;
}


}	// end namespace ilp
}	// end namespace netopt

#endif
