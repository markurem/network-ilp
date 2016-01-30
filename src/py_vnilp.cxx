#include <memory>
#include <vector>
#include <set>

// boost
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/overloads.hpp>

// cplex
#include "andres/ilp/ilpcc.hxx"
#include "andres/ilp/cplex.hxx"

// NetworkILP
#include "NetworkILP.hxx"


/// wrapper utilities.
namespace netopt {
namespace python {

template<class T>
struct VecToList
{
	static PyObject* convert(const std::vector<T>& vec)
	{
		boost::python::list* l = new boost::python::list();
		for(size_t i = 0; i < vec.size(); i++) {
			(*l).append(vec[i]);
		}
		return l->ptr();
	}
};

template<class T>
struct SetToList
{
	static PyObject* convert(const std::set<T>& set)
	{
		boost::python::list* l = new boost::python::list();
		for(typename std::set<T>::iterator iter = set.begin();
				iter != set.end(); ++iter) {
			(*l).append(*iter);
		}
		return l->ptr();
	}
};


/// conversion from python iterable types.
///
struct IterableConverter
{
	/// register converter from python iterable to the provided type.
	template<typename T>
	IterableConverter&
	from_python() {
		boost::python::converter::registry::push_back(
				&IterableConverter::convertible,
				&IterableConverter::construct<T>,
				boost::python::type_id<T>());
		return *this;
	}

	/// check if PyObject is convertible.
	static void* convertible(PyObject* object) {
		return PyObject_GetIter(object) ? object : NULL;
	}

	/// actual converter of PyObject to C++ container
	///
	/// Container concept requirements:
	/// 	* Container::value_type is CopyConstructable
	/// 	* Container can be constructed and populated with two iterators.
	///
	template<typename T>
	static void construct(
			PyObject* object ,
			boost::python::converter::rvalue_from_python_stage1_data* data
	) {
		boost::python::handle<> handle(boost::python::borrowed(object));
		typedef boost::python::converter::rvalue_from_python_storage<T> StorageType;
		void * storage = reinterpret_cast<StorageType*>(data)->storage.bytes;
		typedef boost::python::stl_input_iterator<typename T::value_type> iterType;
		new (storage) T(
				iterType(boost::python::object(handle)),
				iterType());
		data->convertible = storage;
	}
};


}  // end python
}  // end netopt


/// LIBNETOPT
BOOST_PYTHON_MODULE(libnetopt) {
	using namespace boost::python;

	{
		typedef andres::ilp::Cplex<double> ILP;
		/// General Network ILP using CPLEX and CC partitioning.
		{
			typedef andres::ilp::ILPCC<ILP> ILPCC;
			typedef netopt::ilp::NetworkILP<ILPCC> NILPCC;

			// helper for overloaded functions
			typename ILP::Value (NILPCC::*optimizeX1)() = &NILPCC::optimize;
			typename ILP::Value (NILPCC::*optimizeX2)(const bool) = &NILPCC::optimize;

			// class wrapper.
			class_<NILPCC>("NetworkILPCC", init<>())
			.def(init<size_t>())
			.def("add_node", &NILPCC::addNode)
			.def("add_edge", &NILPCC::addEdge)
			.def("add_nodes", &NILPCC::addNodes)
			.def("add_edges", &NILPCC::addEdges)
			.def("add_directed_edges", &NILPCC::addDirectedEdges)
			.def("add_pairwise_weight", &NILPCC::addPairwiseWeight)
			.def("add_pairwise_weights", &NILPCC::addPairwiseWeights)
			.def("add_triplet_weight", &NILPCC::addTripletWeight)
			.def("add_triplet_weights", &NILPCC::addTripletWeights)

			.def("update", &NILPCC::update)
			.def("activate_lazy_connectivity",
				&NILPCC::activateLazyConnectivityConstraints)
			.add_property("root", &NILPCC::getRoot)
			.def("set_root", &NILPCC::setRoot)
			.def("optimize", optimizeX1)
			.def("optimize", optimizeX2)
			.add_property("values", &NILPCC::getValues)
			.add_property("_all_values", &NILPCC::getAllValues)
			.add_property("absolute_gap", &NILPCC::getAbsoluteGap,
				&NILPCC::setAbsoluteGap)
			.add_property("relative_gap", &NILPCC::getRelativeGap,
				&NILPCC::setRelativeGap)
			.add_property("number_of_nodes", &NILPCC::getNumberOfNodes)
			.add_property("number_of_edges", &NILPCC::getNumberOfEdges)
			.add_property("number_of_variables", &NILPCC::getNumberOfVariables)
			.add_property("number_of_constraints", &NILPCC::getNumberOfConstraints)
			;
		}

		/// converters
		to_python_converter<std::vector<typename ILP::Value,
		class std::allocator<typename ILP::Value> >,
		netopt::python::VecToList<typename ILP::Value> >();

		to_python_converter<std::vector<size_t,
		class std::allocator<size_t> >,
		netopt::python::VecToList<size_t> >();

		to_python_converter<std::set<size_t>,
		netopt::python::SetToList<size_t> >();

		netopt::python::IterableConverter()
		/// conversion from iterable.
		.from_python<std::vector<typename ILP::Value> >()
		.from_python<std::vector<float> >()
		.from_python<std::vector<size_t> >()
		/// conversion from nested iterables.
		.from_python<std::vector<std::vector<typename ILP::Value> > >()
		.from_python<std::vector<std::vector<size_t> > >()
		;
	}
}
