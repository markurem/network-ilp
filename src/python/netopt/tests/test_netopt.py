import unittest
from netopt import NetworkILPCC
import numpy as np


class Test(unittest.TestCase):

    def test_add_nodes(self):
        '''addition of decision entities to the model.
        '''
        ilp = NetworkILPCC()
        ilp.add_node(1.3)
        ilp.add_node(-1.7812)
        ilp.add_node(0)
        ilp.add_node(-10)

        assert ilp.number_of_nodes == 4
        assert ilp.number_of_edges == 0

        ilp.add_nodes(np.random.randn(10))
        assert ilp.number_of_nodes == 14
        assert ilp.number_of_edges == 0

    def test_add_edges(self):
        '''insertion of edges to represent neighbouring
        decision entities to the model.
        '''

        ilp = NetworkILPCC()
        ilp.add_nodes(np.random.randn(10))

        ilp.add_edge(0, 1)
        assert ilp.number_of_edges == 1
        ilp.add_edge(1, 0)
        assert ilp.number_of_edges == 2
        ilp.add_edges([3, 1, 2], [1, 2, 3])
        assert ilp.number_of_edges == 8

    def test_basic_setup(self):
        '''basic setup of a model.
        '''
        ilp = NetworkILPCC(3)
        ilp.add_node(5.0)
        ilp.add_node(-2.13)
        ilp.add_node(4)

        ilp.add_nodes([-1, 1, 13.5, 7])

        assert ilp.number_of_nodes == 7
        assert ilp.number_of_edges == 0

        ilp.add_edges([0, 1, 2], [1, 2, 0])
        assert ilp.number_of_edges == 6

        ilp.update()
        assert ilp.number_of_variables == 7

        ilp.add_node(-1)
        ilp.update()
        ilp.number_of_variables == 8
        ilp.number_of_nodes == 8

        assert ilp.root == 0
        ilp.set_root(5)
        assert ilp.root == 5

        val = ilp.optimize(False)
        sol = ilp.values
        assert sol == [0, 1, 0, 1, 0, 0, 0, 1]

        ilp.add_directed_edges([1, 2], [6, 6])
        assert ilp.number_of_edges == 6 + 2

    def test_connectivity(self):
        '''enforced connectivity with constraint generation. 
        '''
        # construct problem.
        node_weights = np.array([[-1, -1, 1, 1],
                                 [-1, -1, 0.5, -1],
                                 [1, 1, 1, -1],
                                 [1, -1, -1, -1]])
        neighbours = np.array([[0,  1], [0,  4], [1,  2], [1,  5], [2,  3],
                               [2,  6], [3,  7], [4,  5], [4,  8], [5,  6],
                               [5,  9], [6,  7], [6, 10], [7, 11], [8,  9],
                               [8, 12], [9, 10], [9, 13], [10, 11], [10, 14],
                               [11, 15], [12, 13], [13, 14], [14, 15]])
        expected = np.array([[1, 1, 0, 0],
                             [1, 1, 1, 1],
                             [0, 0, 0, 1],
                             [0, 1, 1, 1]])
        
        ilp = NetworkILPCC()
        ilp.add_nodes(node_weights.flatten())
        ilp.add_edges(neighbours[:, 0], neighbours[:, 1])
        ilp.update()
        ilp.activate_lazy_connectivity()
        ilp.optimize()
        sol = np.asarray(ilp.values).reshape(node_weights.shape)
        assert np.all(sol == expected)

    def test_higher_order(self):
        '''addition of second and third order weights.
        '''
        ilp = NetworkILPCC()
        ilp.add_nodes([-1, -1, 0.5, -1, 1, -0.5])
        ilp.add_edges([0, 1, 1, 2, 4, 2, ],
                      [1, 2, 4, 4, 5, 3, ], )
        ilp.update()
        ilp.optimize()
        assert np.all(ilp.values == [1, 1, 0, 1, 0, 1])

        ilp.activate_lazy_connectivity()
        ilp.optimize()
        assert np.all(ilp.values == [1, 1, 1, 1, 0, 0])

        before = ilp.number_of_constraints
        ilp.add_pairwise_weight(0, 1, -0.5)
        ilp.optimize()
        assert before == ilp.number_of_constraints

        ilp.add_pairwise_weights([1, 2, 1, 4],
                                 [2, 3, 4, 5],
                                 [0.5, 0.1, -0.5, -0.5])
        ilp.optimize()
        assert np.all(ilp.values == [1, 1, 0, 0, 1, 1])

        print ilp
        ilp.add_triplet_weight(1, 2, 4, -1)
        ilp.optimize()
        assert np.all(ilp.values == [1, 1, 1, 1, 1, 1])


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
