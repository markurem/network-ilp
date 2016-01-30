'''NetworkILP interface.
'''

from .lib.libnetopt import NetworkILPCC as cNetworkILPCC


class NetworkILPCC(cNetworkILPCC):

    def optimize(self, verbose=False):
        '''optimize current problem.
        '''
        if verbose:
            super(NetworkILPCC, self).optimize(verbose)
        else:
            super(NetworkILPCC, self).optimize()

    def __str__(self):
        '''to_string method for NetworkILP.
        '''
        s = '{}'.format(self.__class__.__name__)
        s += '\n   {:<6} : nodes={}, edges={}'.format(
            'graph', self.number_of_nodes, self.number_of_edges)
        s += '\n   {:<6} : variables={}, constraints={}, '.format(
            'ilp',
            self.number_of_variables,
            self.number_of_constraints)
        s += 'gap=(abs: {:2}, rel: {:1.2e})'.format(
            self.absolute_gap,
            self.relative_gap)
        return s
