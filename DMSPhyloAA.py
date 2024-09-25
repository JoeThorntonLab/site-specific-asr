#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 15:56:04 2021

@author: yeonwoo
"""

import getopt, sys
import os
import copy
import numpy as np
from scipy import stats, linalg, optimize
import dendropy


class LogArray:
    '''Log-scale representation of `numpy.ndarray` numeric objects.

    Represents a real number as (log-absolute value, sign) pair and supports
    basic operations that avoid floating point under- or overflow. Negative
    numbers and zero are supported, but infinities are not supported.

    The following operators have been overloaded with support for numpy's
    broadcasting rules:

        `[]`: Indexing and assignment.
        `+`: Element-wise summation.
        `-`: Element-wise subtraction. As a unary operator, negation.
        `*`: Element-wise multiplication.
        `/`: Element-wise division.
        `**`: Element-wise power operation.
        `@`: Matrix (and vector) multiplication.

    Additional operations are supported by the following methods:

        `convert()`
            Back transformation to a regular `numpy.ndarray` object.

        `sum(axis=None, w=None)`
             Weighted sum of array elements over a given axis. Fully
             implemented only for vectors and matrices. For 3-dimensional
             arrays, weighting is not supported. Higher dimensional arrays are
             not supported at all.

        `prod(axis=None)`
            Product of array elements over a given axis. Fully implemented for
            all dimensions.

        `mean(axis=None)`
            Mean of array elements over a given axis. Implemented only for
            3- or lower-dimensional arrays.

        `transpose()`
            Transpose.

        `exp(base=None)`
            Exponentiation.

        `outer(operation='*')`
            Matrix outer-product. Can also perform outer-sum, subtraction,
            and division by specifying operation = '+', '-', or '/'.

    See `__init__` for how to create a LogArray object.
    Attributes
        `logabs` (`numpy.ndarray`)
            Natural log of absolute values of real numbers.
        `sign` (`numpy.ndarray`)
            Signs of real numbers.
        `maxlogabs` (`float`)
            Maximum of `logabs`.
        `ndim` (`int`)
            Dimension of `logabs` and `sign`.
        `shape` (`tuple`)
            Shape of `logabs` and `sign`.
    '''

    __array_ufunc__ = None

    def __init__(self, value, sign=None, isLog=False):
        '''Create and initialize a LogArray object.

        Parameters
        ----------
        value : numpy.ndarray-like
            If isLog is True or sign is given, value must be in log scale.
        sign : numpy.ndarray-like, optional
            Value signs. When provided, value is assumed to be in log scale.
            If sign is None and isLog is True, sign is assumed to be positive.
            The default is None.
        isLog: bool, optional
            Indicates whether value is in log scale. The default is False.
        '''
        if isLog or sign is not None:
            if not isinstance(value, np.ndarray):
                self.logabs = np.array(value)
            else:
                self.logabs = value.copy()
            self.sign = np.ones(self.logabs.shape)
            if sign is not None:
                if self.sign.ndim > 0:
                    self.sign[:] = np.sign(sign)
                else:
                    self.sign = np.sign(sign)
        else:
            self.sign = np.sign(value)
            with np.errstate(divide = 'ignore'):
                self.logabs = np.log(np.abs(value))
        self.maxlogabs = np.max(self.logabs)
        self.ndim = self.logabs.ndim
        self.shape = self.logabs.shape

    def __repr__(self):
        return self.logabs.__repr__() + '\n' + self.sign.__repr__()

    def __getitem__(self, key):
        return LogArray(self.logabs[key], self.sign[key])

    def __setitem__(self, key, other):
        if not isinstance(other, LogArray):
            other = LogArray(other)
        self.logabs[key] = other.logabs.copy()
        self.sign[key] = other.sign.copy()
        self.maxlogabs = np.max(self.logabs)

    def __neg__(self):
        return LogArray(self.logabs.copy(), -self.sign)

    def __abs__(self):
        return LogArray(self.logabs.copy(), np.abs(self.sign))

    def __add__(self, other):
        if not isinstance(other, LogArray):
            other = LogArray(other)
        element_wise_max = np.maximum(self.logabs, other.logabs)
        with np.errstate(invalid = 'ignore', divide = 'ignore'):
            s = ((np.exp(self.logabs - element_wise_max) * self.sign) +
                 (np.exp(other.logabs - element_wise_max) * other.sign))
            if np.isnan(np.sum(s)):
                s[np.isnan(s)] = 0.
            return LogArray(np.log(np.abs(s)) + element_wise_max, np.sign(s))

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if not isinstance(other, LogArray):
            other = LogArray(other)
        element_wise_max = np.maximum(self.logabs, other.logabs)
        with np.errstate(invalid = 'ignore', divide = 'ignore'):
            s = ((np.exp(self.logabs - element_wise_max) * self.sign) -
                 (np.exp(other.logabs - element_wise_max) * other.sign))
            if np.isnan(np.sum(s)):
                s[np.isnan(s)] = 0.
            return LogArray(np.log(np.abs(s)) + element_wise_max, np.sign(s))

    def __rsub__(self, other):
        return self.__sub__(other) * -1.

    def __mul__(self, other):
        if not isinstance(other, LogArray):
            other = LogArray(other)
        return LogArray(self.logabs + other.logabs, self.sign * other.sign)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, LogArray):
            other = LogArray(other)
        return LogArray(self.logabs - other.logabs, self.sign * other.sign)

    def __rtruediv__(self, other):
        if not isinstance(other, LogArray):
            other = LogArray(other)
        return other.__truediv__(self)

    def __matmul__(self, other):
        if not isinstance(other, LogArray):
            other = LogArray(other)
        row_max = np.array(np.max(self.logabs, 1 if self.ndim == 2 else 0))
        row_max[row_max == -np.inf] = 0.
        col_max = np.array(np.max(other.logabs, 0))
        col_max[col_max == -np.inf] = 0.
        if self.ndim == 2:
            res = ((np.exp(self.logabs - row_max[:, np.newaxis]) * self.sign) @
                   (np.exp(other.logabs - col_max) * other.sign))
        else:
            res = ((np.exp(self.logabs - row_max) * self.sign) @
                   (np.exp(other.logabs - col_max) * other.sign))
        with np.errstate(divide = 'ignore'):
            return LogArray(np.log(np.abs(res)) + np.add.outer(row_max, col_max),
                            np.sign(res), isLog=True)

    def __rmatmul__(self, other):
        other = LogArray(other)
        return other.__matmul__(self)

    def __pow__(self, value):
        if isinstance(value, int):
            return LogArray(self.logabs * value, self.sign ** value)
        elif isinstance(value, float):
            if np.any(self.sign < 0):
                raise TypeError('Power of real numbers undefined for negative exponents.')
            return LogArray(self.logabs * value, self.sign)
        else:
            raise TypeError('Power operation supported only for an int or float exponent.')

    def convert(self):
        '''Convert to a numpy.ndarray object.'''
        return np.exp(self.logabs) * self.sign

    def copy(self):
        '''Return a deep copy.'''
        return LogArray(self.logabs, self.sign)

    def exp(self, base=None):
        '''Exponentiate.

        Parameters
        ----------
        base : numeric, optional
            Base to use for exponentiation. Default is e.
        '''
        if base is None:
            base = np.exp(1.)
        return LogArray(self.convert() * np.log(base),
                        sign=np.ones(self.shape), isLog=True)

    def sum(self, axis=None, w=None):
        '''Calculate the sum, optionally along a given axis or weighted.

        Parameters
        ----------
        axis: int, optional
            Axis to be collapsed. The default is None (sum all elements).
        w : numpy.ndarray-like, optional
            Weights to be applied. numpy's broadcasting rule is applied.
            The default is None (= 1.)
        '''
        # for numbers, vectors, or total sums
        if axis is None or self.ndim < 2:
            if self.maxlogabs is np.NINF:
                return LogArray(0.)
            if w is None:
                s = np.sum(np.exp(self.logabs - self.maxlogabs) * self.sign)
            else:
                s = np.sum(np.exp(self.logabs - self.maxlogabs) * self.sign * w)
            return LogArray(np.log(np.abs(s)) + self.maxlogabs, np.sign(s))

        # sliced sums for matrices
        if self.ndim == 2:
            max_vector = np.max(self.logabs, axis)
            with np.errstate(invalid = 'ignore', divide = 'ignore'):
                if axis == 0:
                    if w is None:
                        res = np.sum(np.exp(self.logabs - max_vector) * self.sign, 0)
                    else:
                        res = np.sum(np.exp(self.logabs - max_vector) * self.sign * w[:, np.newaxis], 0)
                elif axis == 1:
                    if w is None:
                        res = np.sum(np.exp(self.logabs - max_vector[:, np.newaxis]) * self.sign, 1)
                    else:
                        res = np.sum(np.exp(self.logabs - max_vector[:, np.newaxis])  * self.sign * w, 1)
                if np.isnan(np.sum(res)):
                    res[np.isnan(res)] = 0.
                return LogArray(np.log(np.abs(res)) + max_vector, np.sign(res))

        # sliced sums for 3D arrays
        if self.ndim == 3:
            if w is not None:
                raise TypeError('Weights are not supported for summing slices of 3D arrays.')
            max_matrix = np.max(self.logabs, axis)
            with np.errstate(invalid = 'ignore', divide = 'ignore'):
                if axis == 0:
                    res = np.sum(np.exp(self.logabs - max_matrix[np.newaxis, ]) * self.sign, 0)
                elif axis == 1:
                    max_matrix_expanded = max_matrix.reshape(max_matrix.shape[0], 1, max_matrix.shape[1])
                    res = np.sum(np.exp(self.logabs - max_matrix_expanded) * self.sign, 1)
                elif axis == 2:
                    res = np.sum(np.exp(self.logabs - max_matrix[:, :, np.newaxis]) * self.sign, 2)
                if np.isnan(np.sum(res)):
                    res[np.isnan(res)] = 0.
                return LogArray(np.log(np.abs(res)) + max_matrix, np.sign(res))

        raise TypeError('Method .sum() is not supported for >3D arrays.')

    def prod(self, axis=None):
        '''Calculate the product, optionally along a given axis.'''
        return LogArray(np.sum(self.logabs, axis), np.prod(self.sign, axis))

    def mean(self, axis=None):
        '''Calculate the mean, optionally along a given axis.'''
        if axis is None or self.ndim < 2:
            if self.maxlogabs is np.NINF:
                return LogArray(0.)
            else:
                res = np.mean(np.exp(self.logabs - self.maxlogabs) * self.sign)
                return LogArray(np.log(np.abs(res)) + self.maxlogabs, np.sign(res))

        # sliced mean for matrices
        if self.ndim == 2:
            max_vector = np.max(self.logabs, axis)
            with np.errstate(invalid = 'ignore', divide = 'ignore'):
                if axis == 0:
                    res = np.mean(np.exp(self.logabs - max_vector) * self.sign, 0)
                elif axis == 1:
                    res = np.mean(np.exp(self.logabs - max_vector[:, np.newaxis]) * self.sign, 1)
                if np.isnan(np.sum(res)):
                    res[np.isnan(res)] = 0.
                return LogArray(np.log(np.abs(res)) + max_vector, np.sign(res))

        # sliced mean for 3D arrays
        if self.ndim == 3:
            max_matrix = np.max(self.logabs, axis)
            with np.errstate(invalid = 'ignore', divide = 'ignore'):
                if axis == 0:
                    res = np.mean(np.exp(self.logabs - max_matrix[np.newaxis, ]) * self.sign, 0)
                elif axis == 1:
                    max_matrix_expanded = max_matrix.reshape(max_matrix.shape[0], 1, max_matrix.shape[1])
                    res = np.mean(np.exp(self.logabs - max_matrix_expanded) * self.sign, 1)
                elif axis == 2:
                    res = np.mean(np.exp(self.logabs - max_matrix[:, :, np.newaxis]) * self.sign, 2)
                if np.isnan(np.sum(res)):
                    res[np.isnan(res)] = 0.
                return LogArray(np.log(np.abs(res)) + max_matrix, np.sign(res))

        raise TypeError('Method .mean() is not supported for >3D arrays.')

    def transpose(self):
        '''Return the transpose.'''
        if self.ndim < 2:
            return self.copy()
        elif self.ndim == 2:
            return LogArray(self.logabs.T, self.sign.T)
        else:
            raise TypeError('Transpose is supported only for numbers, vectors, and matrices.')

    def outer(self, operation='*'):
        '''Generate an outer product of a vector,

        Parameters
        ----------
        operation : string, optional
            Operation to be performed on each element. The default is
            multiplication ('*'). Supports addition ('+'), subtraction ('-'),
            and division ('/').
        '''
        if(self.ndim > 1):
            raise TypeError('Outer operation only supported for vectors.')
        tiled = LogArray(np.tile(self.logabs, (self.shape[0], 1)),
                         sign=np.tile(self.sign, (self.shape[0], 1)),
                         isLog=True)
        if operation == '*':
            return tiled.transpose() * tiled
        elif operation == '+':
            return tiled.transpose() + tiled
        elif operation == '-':
            return tiled.transpose() - tiled
        elif operation == '/':
            return tiled.transpose() / tiled
        else:
            raise TypeError('Invalid operation.')

class DMSPhyloAA:
    '''Deep-mutational scanning-informed protein phylogenetics.

    See `__init__` for how to initialize a `DMSPhyloAA` object.


    Example usage

    1. Creating and initializing a DMSPhyloAA object.
    model = DMSPhyloAA(alignmentFile, treeFile, DMSFile,
                       DNAModel='GTR', fitnessFuncType='logistic', ASRVk=4,
                       DNAR=None, DNAPI=None, fitnessFuncParam=None, ASRVa=1)

    2. Optimizing and exporting independent parameters.
    model.optimize()
    model.export(dir)

    3. Performing ASR for all nodes and exporting the results.
    model.export_ASR(dir)

    4. Ignoring DMS data and using the uniform model of protein evolution.
    model._derive_protein_model(R=np.ones((model.nsite, 190)),
                                PI=np.ones((model.nsite, 20)))
    model._calc_llk()
    model.optimize(optProteinModel=False)


    Attributes

    Data:
        `alignment` (`numpy.ndarray` of shape `(ntaxa, nsite)`)
            Sequence alignment. Amino acids are coded as integers from 0 to 19
            ordered as in the attribute `AA`. Gaps are coded as -1.
        `tree` (instance of `dendropy.Tree`)
            Phylogenetic tree. Tree nodes have attributes `llk`, the
            conditional log-likelihood profile, and 'label', a node label.
        `DMS` (`numpy.ndarray` of shape `(nsite, 20)`)
            Mutant phenotypes. No missing data are allowed.

    Model specifiers:
        `DNAModel` (`str`)
            DNA model. Currently one of 'JC69', 'HKY85' and 'GTR'.
        `fitnessFuncType` (`str`)
            Parametric function converting mutant phenotypes into growth rates
            (called 'fitness function'). Currently one of 'logistic'.
        `ASRVk` (`int`)
            Number of categories of discrete gamma distribution for modeling
            among-site rate variation (ASRV). Set to 1 to disable ASRV.

    Independent model parameters:
        `DNAR` (`numpy.ndarray` of shape `(4, 4)`)
            Nucleotide exchangeability matrix.
        `DNAPI` (`numpy.ndarray` of shape `(4,)`)
            Nucleotide equilibrium frequency vector.
        `fitnessFuncParam` (`numpy.ndarray`)
            Fitness function parameter vector; the length and values depend on
            `fitnessFuncType`.
        `ASRVa` (`float`)
            Discrete gamma distribution shape parameter.

    Dependent model parameters:
        `codonPI` (`numpy.ndarray` of shape `(64,)`)
            Equilibrium codon frequencies implied by the mutation process.
        `MU` (`numpy.ndarray` of shape `(20, 20)`)
            Amino acid mutation rate matrix.
        `S` (`numpy.ndarray` of shape `(nsite, 20, 20)`)
            Site-specific amino acid selection coefficient matrices.
        `R` (`numpy.ndarray` of shape `(nsite, 20, 20)`)
            Site-specific amino acid exchangeability matrices.
        `PI` (`numpy.ndarray` of shape `(nsite, 20)`)
            Site-specific amino acid equilibrium frequency vectors.
        `PINoSelection` (`numpy.ndarray` of shape `(20,)`)
            Amino acid equilibrium frequencies implied by the mutation process.
        `Q` (`numpy.ndarray` of shape `(nsite, 20, 20)`)
            Site-specific amino acid transition rate matrices.
        `LQ` (`numpy.ndarray` of shape `(nsite, 20, 20)`)
            Left eigenvector matrices of Q
        `RQ` (`numpy.ndarray` of shape `(nsite, 20, 20)`)
            Right eigenvector matrices of Q
        `AQ` (`numpy.ndarray` of shape `(nsite, 20)`)
            Eigenvalue vectors of Q
        `relRate` (`numpy.ndarray` of shape `(nsite,)`)
            Relative rate of evolution across sites.
        `ASRVr` (`numpy.ndarray` of shape `(ASRVk,)`)
            Relative rate for each discrete gamma distribution category.

    Likelihood-related:
        `llkSite` (`numpy.ndarray` of shape `(nsite,)`)
            Log-likelihood for each site.
        `llk` (`numpy.ndarray` of shape `(1,)`)
            Total log-likelihood.

    Constants:
        `nsite` (`int`)
            Number of sites in the alignment.
        `ntaxa` (`int`)
            Number of taxa in the alignment and tree.
        `taxa` (`list` of `str`)
            Taxon names. Must be identical between the alignment and tree.
        `AA` (`list` of `str`)
            Amino acid symbols. The order is identical to that in most
            phylogenetics softwares.
        `DNA` (`list` of `str`)
            Nucleotide symbols.
        `GENETIC_CODE` (`tuple`)
            Standard genetic code.

    Control parameters:
        `_minbl` (`float`)
            Lower bound of branch length.
        '_maxbl' (`float`)
            Upper bound of branch length.
        `_mingr` (`float`)
            Lower bound of growth rate.
        '_maxbl' (`float`)
            Upper bound of growth rate.
        `_minval` (`float`)
            L-BFGS-B lower bound for all other parameters that must be positive.
        `_tol` (`float`)
            Total log-likelihoood accuracy goal for optimization.
    '''

    AA = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I','L', 'K', 'M', 'F',
          'P', 'S', 'T', 'W', 'Y', 'V')
    DNA = ('A', 'C', 'G', 'T')
    GENETIC_CODE = (('AAA', 'K'),('AAC', 'N'),('AAG', 'K'),('AAT', 'N'),
                    ('ACA', 'T'),('ACC', 'T'),('ACG', 'T'),('ACT', 'T'),
                    ('AGA', 'R'),('AGC', 'S'),('AGG', 'R'),('AGT', 'S'),
                    ('ATA', 'I'),('ATC', 'I'),('ATG', 'M'),('ATT', 'I'),
                    ('CAA', 'Q'),('CAC', 'H'),('CAG', 'Q'),('CAT', 'H'),
                    ('CCA', 'P'),('CCC', 'P'),('CCG', 'P'),('CCT', 'P'),
                    ('CGA', 'R'),('CGC', 'R'),('CGG', 'R'),('CGT', 'R'),
                    ('CTA', 'L'),('CTC', 'L'),('CTG', 'L'),('CTT', 'L'),
                    ('GAA', 'E'),('GAC', 'D'),('GAG', 'E'),('GAT', 'D'),
                    ('GCA', 'A'),('GCC', 'A'),('GCG', 'A'),('GCT', 'A'),
                    ('GGA', 'G'),('GGC', 'G'),('GGG', 'G'),('GGT', 'G'),
                    ('GTA', 'V'),('GTC', 'V'),('GTG', 'V'),('GTT', 'V'),
                    ('TAA', '*'),('TAC', 'Y'),('TAG', '*'),('TAT', 'Y'),
                    ('TCA', 'S'),('TCC', 'S'),('TCG', 'S'),('TCT', 'S'),
                    ('TGA', '*'),('TGC', 'C'),('TGG', 'W'),('TGT', 'C'),
                    ('TTA', 'L'),('TTC', 'F'),('TTG', 'L'),('TTT', 'F'))
    _minbl = 1e-4
    _maxbl = 10.
    _mingr = 1e-3
    _maxgr = 15.
    _minval = 1e-3
    _tol = 1e-4

    def __init__(self, alignmentFile, treeFile, DMSFile, DNAModel,
                 fitnessFuncType, ASRVk, DNAR=None, DNAPI=None,
                 fitnessFuncParam=None, ASRVa=1):
        '''Create and initialize a DMSPhyloAA object.

        Parameters
        ----------
        alignmentFile : str
            Path to a phylip-format alignment file.
        treeFile : str
            Path to a newick-format tree file.
        DMSFile : str
            Path to a csv-format DMS data file.
        DNAModel : str
            DNA model type; must be one of 'JC69', 'HKY85', or 'GTR'.
        DNAR : numeric(6), optional
            Lower triangular part of the initial DNA R matrix.
            The default is all 1.
        DNAPI : numeric(4), optional
            Initial equilibrium base frequencies. The default is all 0.25.
        fitnessFuncType : str
            Fitness function type; must be one of 'logistic'.
        fitnessFuncParam : numeric, optional
            Initial fitness function parameters.
        ASRVk : int
            Number of discrete gamma distribution categories for modeling
            among-site rate variation (ASRV). Use 1 to disable ASRV.
        ASRVa : numeric(1), optional
            Initial ASRV shape parameter. The default is 1.
        '''
        self._import_alignment(alignmentFile)
        self._import_tree(treeFile)
        self._import_DMS(DMSFile)
        self._set_DNAModel(DNAModel, DNAR, DNAPI)
        self._set_fitness_function(fitnessFuncType, fitnessFuncParam)
        self._set_ASRV(ASRVk, ASRVa)

        # Calculating the conditional likelihood profiles of terminal nodes
        for node in self.tree.postorder_node_iter():
            node.llk = np.zeros((self.nsite, self.ASRVk, 20))
            if node.is_leaf():
                for i in range(self.nsite):
                    aa = self.alignment[self.taxa.index(node.taxon.label), i]
                    if aa != -1:
                        node.llk[i,] = np.full((self.ASRVk, 20), np.NINF)
                        node.llk[i, :, aa] = 0.

        self._derive_protein_model()
        self._calc_llk()

    def optimize(self, optProteinModel=True, optBl=True, optASRV=True,
                 useInitBl=False, maxIter=30):
        '''Optimize the independent model parameters.

        Parameters
        ----------
        optProteinModel : bool, optional
            Should the DNA model and fitness function be optimized?
            The default is True.
        optBl : bool, optional
            Should branch lengths be optimized? The default is True.
        optASRV : bool, optional
            Should the ASRV shape parameter be optimized? The default is True.
        nrounds : int, optional
            Maximum number of rounds of optimization. The default is 20.
        useInitBl : bool, optional
            Should branch lengths in the tree file be used as initial values?
            If False, all branch lengths are set to the average of branch
            lengths in the tree file. The default is False.
        '''
        print('Optimization')
        print(f'DNA model: {self.DNAModel}')
        print(f'Fitness function: {self.fitnessFuncType}')
        print(f'Number of ASRV categories: {self.ASRVk}')

        # A round of crude optimization
        llkInit = self.llk
        print(f'Round 1, crude optimization: initial llk = {llkInit:.6f}.')
        if optProteinModel:
            self._optimize_protein_model(1e-3)
        if optBl:
            if not useInitBl:
                bl = np.delete([edge.length for edge in self.tree.edges()], 0)
                avgbl = np.mean(bl)
                self._set_uniform_bl(avgbl)
                self._calc_llk()
                print(f'Setting all initial branch lengths to {avgbl:.4f}; ' +
                      f'llk = {self.llk:.6f}')
            self._optimize_bl(1.)
        if optASRV:
            self._optimize_ASRV(1e-3)
        print(f'End of round 1: initial llk = {llkInit:.6f}, ' +
              f'final llk = {self.llk:.6f}.')

        for i in range(1, maxIter):
            llkInit = self.llk
            print(f'Round {i+1}: initial llk = {llkInit:.6f}.')
            if optProteinModel:
                self._optimize_protein_model()
            if optBl:
                self._optimize_bl()
            if optASRV:
                self._optimize_ASRV()
            print(f'End of round {i+1}: initial llk = {llkInit:.6f}, ' +
                  f'final llk = {self.llk:.6f}')
            if self.llk < llkInit + self._tol:
                print('Optimization complete!')
                return None
        print('Maximum number of iteration reached.')
        print('Rerun using useInitBl=True for further optimization.')

    def ASR(self, node):
        '''Ancestral sequence reconstruction for a node. Currently only
        implemented for bifurcating trees.

        Parameters
        ----------
        node : int
            Node label. Check the exported labeled cladogram for node labels.

        Returns
        -------
        MPA : str
            Maximum a posteriori sequence.
        PP : np.ndarray (nsite, 20)
            Posterior probability distribution.
        '''
        node = self.tree.find_node_with_label(node)
        if node.is_leaf():
            raise ValueError(f'Node {node.label} is a terminal node.')
        if node == self.tree.seed_node and self.tree.root_state == 2:
            raise ValueError(f'Node {node.label} is the root node.')

        # Conditional likelihood profiles of adjacent subtrees
        if node == self.tree.seed_node:
            L1, L2, L3 = (LogArray(child.llk, isLog=True) for child in node.child_node_iter())
            b1, b2, b3 = (child.edge_length for child in node.child_node_iter())
        else:
            L1, L2 = (LogArray(child.llk, isLog=True) for child in node.child_node_iter())
            b1, b2 = (child.edge_length for child in node.child_node_iter())
            L3 = LogArray(self._calc_subtree_llk(node), isLog=True)
            if node.parent_node == self.tree.seed_node and self.tree.root_state == 2:
                b3 = 2. * node.edge_length
            else:
                b3 = node.edge_length

        PP = np.zeros((self.nsite, 20))
        Liks = LogArray(np.ones((self.nsite, self.ASRVk, 20)))
        for i in range(self.nsite):
            for k in range(self.ASRVk):
                P1 = self._calc_P(i, b1, k)
                P2 = self._calc_P(i, b2, k)
                P3 = self._calc_P(i, b3, k)
                Liks[i, k, ] = ((P1 @ L1[i, k, ]) * (P2 @ L2[i, k, ]) *
                                (P3 @ L3[i, k, ]))
            rawPP = Liks[i,].mean(axis=0) * self.PI[i,]
            PP[i,] = (rawPP / rawPP.sum()).convert()

        MPA = ''
        for i in np.argmax(PP, axis=1):
            MPA += self.AA[i]
        return MPA, PP

    def update_model(self):
        '''Recalculate the dependent model parameters and likelihood using the
        current values of independent parameters.'''
        self._set_DNAModel(self.DNAModel, self.DNAR[np.tril_indices(4, -1)],
                           self.DNAPI)
        self._set_fitness_function(self.fitnessFuncType, self.fitnessFuncParam)
        self._set_ASRV(self.ASRVk, self.ASRVa)
        self._derive_protein_model()
        return self._calc_llk()

    def export(self, dirPath):
        '''Export the tree and independent model parameters.
        Parameters
        ----------
        dirPath : str
            Directory to which text files should be written. Must end with /.
        '''
        try:
            self.tree.write(path=dirPath+'phylogram.txt', schema='newick',
                            suppress_rooting=True,
                            suppress_internal_node_labels=True)
        except:
            raise ValueError('Provide a valid directory path.')
        self.tree.write(path=dirPath+'labeled_cladogram.txt', schema='newick',
                        suppress_rooting=True, suppress_edge_lengths=True)
        with open(dirPath+'inferred_parameters.txt', 'w') as f:
            f.writelines([f'DNA model: {self.DNAModel}\n',
                          '   R (lower triangular; row-by-row): ',
                          np.array2string(self.DNAR[np.tril_indices(4, -1)]),
                          '\n   PI: ',
                          np.array2string(self.DNAPI),
                          f"\nASRV shape: {self.ASRVa if self.ASRVk > 1 else 'NA'}",
                          f' (number of categories: {self.ASRVk})\n',
                          f'Fitness function ({self.fitnessFuncType}) parameters: ',
                          np.array2string(self.fitnessFuncParam),
                          f'\ntotal log-likelihood: {self.llk:.4f}\n'
                          ])
        np.savetxt(dirPath+'per_site_llk.txt', self.llkSite, fmt='%.4f')
        np.savetxt(dirPath+'AAPI.txt', self.PI, fmt='%.4f')

    def export_ASR(self, dirPath, nodeList=None):
        '''Export MPA sequences and posterior probability distributions for a
        list of internal nodes.

        Parameters
        ----------
        dirPath : str
            Directory to which text files should be written. Must end with /.
        nodeList : numeric, optional
            A list of node labels. The default is all internal nodes.
        '''
        if nodeList is None:
            nodeList = []
            for node in self.tree.preorder_internal_node_iter():
                if node == self.tree.seed_node and self.tree.root_state == 2:
                    continue
                else:
                    nodeList.append(node.label)
        elif isinstance(nodeList, int): nodeList = [nodeList]

        header = ''
        for i, aa in enumerate(self.AA):
            header += '"' + str(aa) + '"'
            if i < 19: header += ','

        for node in nodeList:
            MPA, PP = self.ASR(node)
            with open(dirPath + str(node) + '.txt', 'w') as file:
                file.write(MPA)
            np.savetxt(dirPath + str(node) + '_PP.txt', PP, delimiter=',',
                       header=header, comments='', fmt='%.3f')

    def export_state_tree(self, dirPath):
        '''Export a list of trees with reconstructed ancestral states shown
        as node labels.

        Parameters
        ----------
        dirPath : str
            Directory to which text files should be written. Must end with /.
        '''
        for node in self.tree.preorder_internal_node_iter():
            if node == self.tree.seed_node and self.tree.root_state == 2:
                continue
            else:
                MPA, PP = self.ASR(node.label)
                node.ASR = []
                for i in range(self.nsite):
                    maxpp = np.max(PP[i,])
                    node.ASR.append(MPA[i] + ' (' + str(maxpp)[0:5] + ')')

        with open(dirPath + 'state_labeled_trees.txt', 'a') as file:
            for site in range(self.nsite):
                tree_to_print = copy.deepcopy(self.tree)
                for node in tree_to_print.preorder_internal_node_iter():
                    if node == tree_to_print.seed_node and tree_to_print.root_state == 2:
                        continue
                    node.label = node.ASR[site]
                for node in tree_to_print.leaf_node_iter():
                    state = self.alignment[self.taxa.index(node.taxon.label), site]
                    if state == -1:
                        state = '-'
                    else:
                        state = self.AA[state]
                    node.taxon.label += ' (' + state + ')'
                tree_to_print.write(file=file, schema='newick', suppress_rooting=True)

    def _import_alignment(self, alignmentFile):
        # Import a phylip format alignment file.
        try:
            with open(alignmentFile) as file:
                self.ntaxa, self.nsite = map(int, file.readline().split())
                self.taxa = [''] * self.ntaxa
                alignmentStr = [''] * self.ntaxa
                for i, line in enumerate(file.readlines()):
                    self.taxa[i], alignmentStr[i] = line.split()
        except:
            raise ValueError('Could not import alignment. ' +
                             'Provide a phylip-format alignment.')

        # Converting amino acid symbols into integers
        self.alignment = np.zeros((self.ntaxa, self.nsite), int)
        for i, seq in enumerate(alignmentStr):
            seqInt = []
            for aa in list(seq):
                if aa in self.AA:
                    seqInt.append(self.AA.index(aa))
                else:
                    seqInt.append(-1)  # gap or missing data
            self.alignment[i, ] = seqInt

    def _import_tree(self, treeFile, schema='newick'):
        # Import a tree file as a dendropy.Tree object.
        try:
            self.tree = dendropy.Tree.get(path=treeFile, schema=schema,
                                          preserve_underscores=True)
        except:
            raise ValueError('Could not import tree. ' +
                             'Provide a newick-format tree.')

        self.tree.reroot_at_node(self.tree.seed_node)
        self.tree.root_state = len(self.tree.seed_node.child_nodes())
        # If a cladogram, all branch lengths are set to 0.1
        if self.tree.seed_node.child_edges()[0].length is None:
            self._set_uniform_bl(0.1)
        # Node labeling; the order of labeling is crucial; do not alter
        for i, node in enumerate(self.tree.preorder_node_iter()):
            node.label = i + 1

    def _import_DMS(self, DMSFile):
        '''Import DMS data.

        Data must be in csv format with rows corresponding to sites and columns to
        amino acids. Single-letter amino acid symbols must be used as column names.
        No row names should be present.
        '''
        if DMSFile == '':
            self.DMS = np.zeros((self.nsite, 20))
        else:
            try:
                DMS = np.genfromtxt(DMSFile, delimiter=',', names=True)
                self.DMS = np.zeros((self.nsite, 20))
                for i, aa in enumerate(self.AA):
                    self.DMS[:, i] = DMS[aa]
                del DMS
            except:
                raise ValueError('Could not import DMS file.\n'+
                                 'Check ._import_DMS() for format requirements.')

    def _set_DNAModel(self, DNAModel, R, PI):
        # Set the DNA model.
        if DNAModel not in {'JC69', 'HKY85', 'GTR'}:
            raise ValueError('Only JC69, HKY85, and GTR are supported.')
        self.DNAModel = DNAModel
        R = np.ones(6) if R is None else np.array(R)
        if len(R) != 6:
            raise ValueError('Provide 6 lower triangular values for DNA R.')
        self.DNAR = np.full((4, 4), np.nan)
        self.DNAR[np.tril_indices(4, -1)] = R
        self.DNAR[([0, 0, 1, 0, 1, 2], [1, 2, 2, 3, 3, 3])] = R
        PI = np.ones(4) if PI is None else np.array(PI)
        if len(PI) != 4:
            raise ValueError('Provide 4 values for DNA PI.')
        self.DNAPI = PI / np.sum(PI)

    def _set_fitness_function(self, fitnessFuncType, fitnessFuncParam):
        # Set the function for converting mutant phenotypes into growth rates.
        if fitnessFuncType not in {'logistic'}:
            raise ValueError('Only logistic function is supported.')
        self.fitnessFuncType = fitnessFuncType
        if fitnessFuncType == 'logistic':
            if fitnessFuncParam is None:
                self.fitnessFuncParam = np.array([1., 1., 0.])
            else:
                fitnessFuncParam = np.array(fitnessFuncParam)
                if len(fitnessFuncParam) != 3:
                    raise ValueError('Provide three values for logistic ' +
                                     'function parameters: ' +
                                     'Maximum, Steepness, Midpoint')
                if fitnessFuncParam[0] > self._maxgr:
                    print('Maximum growth rate modified to be within bound.')
                    fitnessFuncParam[0] = self._maxgr
                self.fitnessFuncParam = fitnessFuncParam

    def _fitness_function(self, x):
        # Function for converting mutant phenotypes into growth rates.
        if self.fitnessFuncType == 'logistic':
            exponent = np.clip(-self.fitnessFuncParam[1] * (x - self.fitnessFuncParam[2]),
                               -100., 100.)
            return self.fitnessFuncParam[0] / (1. + np.exp(exponent)) + self._mingr

    def _set_ASRV(self, k, a):
        # Set the discrete gamma-distributed ASRV model.
        self.ASRVk = k
        self.ASRVa = None if k == 1 else a
        if k == 1:
            self.ASRVr = np.array([1.])
        else:
            q = stats.gamma.ppf(np.array(range(1, k), float)/k, a=a, scale=1/a)
            self.ASRVr = k * np.diff(np.concatenate(([0.], stats.gamma.cdf(q * a, a=a+1), [1.])))

    def _set_R(self, RLowerTri):
        # Set the amino acid R matrices.
        if RLowerTri.shape != (self.nsite, 190):
            raise ValueError('Overriding R: Incorrect input dimensions.')
        self.R = np.zeros((self.nsite, 20, 20))
        for i in range(self.nsite):
            Ri = np.zeros((20, 20))
            Ri[np.tril_indices(20, -1)] = RLowerTri[i, ]
            self.R[i, ] = Ri + Ri.transpose()
            np.fill_diagonal(self.R[i, ], np.nan)

    def _set_PI(self, PI):
        # Set the amino acid PI vectors.
        if PI.shape != (self.nsite, 20):
            raise ValueError('Overriding PI: Incorrect input dimensions.')
        self.PI = np.zeros((self.nsite, 20))
        for i in range(self.nsite):
            self.PI[i, ] = PI[i, ] / np.sum(PI[i, ])

    def _set_uniform_bl(self, bl):
        # Set all branch lengths to bl.
        for node in self.tree.postorder_node_iter():
            if node != self.tree.seed_node:
                node.edge_length = bl
        if self.tree.root_state == 2:
            for child in self.tree.seed_node.child_node_iter():
                child.edge_length = bl / 2.

    def _derive_protein_model(self, R=None, PI=None):
        '''Derive amino acid R and PI for each site.

        R or PI can be specified, which will override origin-fixation model
        derived from the DNA model and DMS data. In this case, R must be a 2D
        numpy.ndarray of dimension (nsite, 190), each row containing the 190
        lower triangular elements of R for a site row-by-row. PI must be a 2D
        numpy.ndarray of dimension (nsite, 20), each row containing the 20
        equilibrium frequency values for a site.
        '''
        if R is not None:
            self._set_R(R)
            if PI is not None:
                self._set_PI(PI)
                self._calc_decompose_Q()
                return None
        elif PI is not None:
            self._set_PI(PI)

        # Integer encoding of codons and amino acids
        codonInt, aaInt = [], []
        for codon, aa in self.GENETIC_CODE:
            codonInt.append([self.DNA.index(base) for base in codon])
            aaInt.append(self.AA.index(aa) if aa != '*' else 20)
        codonInt, aaInt = np.array(codonInt), np.array(aaInt)

        # Equilibrium frequency of codons and amino acids under no selection
        self.codonPI = np.array([np.prod([self.DNAPI[base] for base in codon]) for codon in codonInt])
        self.PINoSelection = np.array([np.sum(self.codonPI[aaInt == i]) for i in range(20)])
        self.PINoSelection /= np.sum(self.PINoSelection)

        # Calculating the mutation rate matrix MU
        self.MU = np.zeros((20, 20))
        np.fill_diagonal(self.MU, np.nan)
        DNAQ = self.DNAR * self.DNAPI
        for i in range(20):
            for j in range(20):
                if i == j:
                    continue
                codonsForI = np.flatnonzero(aaInt == i)
                codonsForJ = np.flatnonzero(aaInt == j)
                # calculating the rate of change from codon i to codon j
                for codonI in codonsForI:
                    for codonJ in codonsForJ:
                        diffPos = np.flatnonzero(codonInt[codonI, ] != codonInt[codonJ, ])
                        if(len(diffPos) != 1):
                            continue
                        self.MU[i, j] += DNAQ[codonInt[codonI, diffPos], codonInt[codonJ, diffPos]] \
                            * self.codonPI[codonI] / np.sum(self.codonPI[codonsForI])

        # Calculating the selection coefficient matrix S
        #   and the probability of fixation matrix PF
        self.S = np.zeros((self.nsite, 20, 20))
        gr = self._fitness_function(self.DMS)
        for i in range(self.nsite):
            Si = -np.subtract.outer(gr[i, ], gr[i, ])
            Si[np.abs(Si) < 1e-5] = 0.
            self.S[i, ] = Si
        with np.errstate(invalid='ignore'):
            self.PF = self.S / (1. - np.exp(-2. * self.S))
            if np.isnan(np.sum(self.PF)):
                self.PF[np.isnan(self.PF)] = 0.5

        # Calculating the amino acid equilibrium frequency vectors
        if PI is None:
            self.PI = np.exp(2. * gr) * self.PINoSelection
            self.PI /= np.sum(self.PI, axis=1)[:, np.newaxis]

        # Calculating the exchangeability matrices
        if R is None:
            self.R = np.zeros((self.nsite, 20, 20))
            for i in range(self.nsite):
                self.R[i,] = self.MU * self.PF[i,] / self.PI[i,]

        # Calculating and decomposing the transition rate matrices
        self._calc_decompose_Q()

    def _calc_decompose_Q(self):
        '''Calculate the transition rate matrix Q and its eigen-decomposition
        for each site.

        Note that the Q matrix for each site is normalized by the mean rate of
        substitution of all sites. This means that changing the Q matrix for
        one site must be accompanied by re-normalization of all sites.
        '''
        self.Q = np.zeros((self.nsite, 20, 20))
        self.LQ = np.zeros((self.nsite, 20, 20)) # left eigenmatrices of Q
        self.RQ = np.zeros((self.nsite, 20, 20)) # right eigenmatrices of Q
        self.AQ = np.zeros((self.nsite, 20)) # eigenvalues of Q
        self.relRate = np.zeros(self.nsite) # relative rate of substitution

        # Calculating the relative rates of evolution and normalized Q
        for i in range(self.nsite):
            self.Q[i,] = self.R[i,] * self.PI[i,]
            np.fill_diagonal(self.Q[i,], -np.nansum(self.Q[i,], axis=1))
            self.relRate[i] = -np.sum(np.diag(self.Q[i,]) * self.PI[i,])
        self.Q /= np.mean(self.relRate)
        self.relRate /= np.mean(self.relRate)

        # Eigen-decomposition of Q
        for i in range(self.nsite):
            PISqrt = np.diag(np.sqrt(self.PI[i,]))
            PISqrtInv = np.diag(np.sqrt(1. / self.PI[i,]))
            B = PISqrt @ self.Q[i,] @ PISqrtInv
            A, U = linalg.eigh(B)
            self.AQ[i,] = A
            self.LQ[i,] = PISqrtInv @ U
            self.RQ[i,] = U.T @ PISqrt

    def _calc_P(self, i, bl, k):
        '''Calculate the transition probability matrix P.

        Due to numerical cancellation during the multiplication of the eigen-
        vectors and values of Q, elements of P may be negative. To prevent
        downstream errors, elements of P smaller than 1e-10 are converted to
        1e-10.

        Parameters
        ----------
        i : int
            Site (numbering begins at 0).
        bl : numpy.ndarray-like
            Branch length.
        k : int
            Rate category (numbering begins at 0.

        Returns
        -------
        P : numpy.ndarray (20, 20)
            Transition probability matrix.
        '''
        P = self.LQ[i,] @ np.diag(np.exp(self.AQ[i,] * self.ASRVr[k] * bl)) @ self.RQ[i,]
        P[P < 1e-10] = 1e-10
        return P

    def _calc_llk(self):
        '''Calculate the conditional log-likelihood profile of every node and
        the total log-likelihood.

        The conditional likelihood profile of a node is the likelihood of the
        subtree rooted at the node calculated for each ancestral state and
        rate cateogory. It is stored as the .llk attribute of each node.
        '''
        for node in self.tree.postorder_internal_node_iter():
            node.llk = np.zeros((self.nsite, self.ASRVk, 20))
            for child in node.child_node_iter():
                node.llk += self._prune_branch(child.llk, child.edge_length)
        self.llkSite = np.zeros(self.nsite)
        lk = LogArray(self.tree.seed_node.llk, isLog=True)
        for i in range(self.nsite):
            self.llkSite[i] = lk[i,].mean(axis=0).sum(w = self.PI[i,]).logabs
        self.llk = np.sum(self.llkSite)
        return self.llk

    def _calc_posterior_rate(self):
        '''Calculate posterior mean rate of evolution for each site.

        The relative rate of evolution determined by the experimental model is
        combined with the posterior mean rate determined by the among-site
        rate variation model.
        '''
        lk = LogArray(self.tree.seed_node.llk, isLog=True)
        prate = np.zeros(self.nsite)
        for i in range(self.nsite):
            lk_i = lk[i,].sum(axis=1, w = self.PI[i,])
            prate[i] = (lk_i * self.ASRVr / lk_i.sum()).sum().convert()
        return prate * self.relRate

    def _prune_branch(self, llk, bl):
        '''Given the conditional log-likelihood profile of a node and branch
        length to its parent node, calculate the conditional log-likelihood
        profile of the parent node.

        Parameters
        ----------
        llk : numpy.ndarray (nsite, ASRVk, 20)
            Conditional log-likelihood profile of a node.
        bl : numpy.ndarray-like
            Branch length to the parent node.
        '''
        llkParent = np.zeros((self.nsite, self.ASRVk, 20))
        for i in range(self.nsite):
            for k in range(self.ASRVk):
                P = self._calc_P(i, bl, k)
                maxllk = np.max(llk[i, k, ])
                normllk = np.exp(llk[i, k, ] - maxllk)
                llkParent[i, k, ] = np.log(P @ normllk) + maxllk
        return llkParent

    def _calc_subtree_llk(self, node):
        '''Calculate the conditional likelihood profile of a subtree generated
        by removing the subtree rooted at the node and rerooting the remaining
        tree at the ancestor of the node.
        '''
        subtree = copy.deepcopy(self.tree)
        node = subtree.find_node_with_label(node.label)
        subtree.reroot_at_node(node.parent_node)
        subtree.prune_subtree(node)

        for node in subtree.postorder_internal_node_iter():
            childs = {child.label for child in node.child_node_iter()}
            nodeOri = self.tree.find_node_with_label(node.label)
            childsOri = {child.label for child in nodeOri.child_node_iter()}
            if childs == childsOri:
                continue
            node.llk = np.zeros((self.nsite, self.ASRVk, 20))
            for child in node.child_nodes():
                node.llk += self._prune_branch(child.llk, child.edge_length)
        return subtree.seed_node.llk

    def _calc_llk_bl(self, L1, L2, bl):
        '''Given the length of a branch and the conditional likelihood profiles
        at the two flanking nodes, calculate the total likelihood.

        Parameters
        ----------
        L1, L2 : LogArray
            Conditional likelihood profiles at the two flanking nodes.
        bl : numpy.ndarray-like
            Branch length.
        '''
        L = LogArray(np.ones((self.nsite, self.ASRVk)))
        for i in range(self.nsite):
            for k in range(self.ASRVk):
                L[i, k] = (self.PI[i,] * L1[i, k, ]) @ (self._calc_P(i, bl, k) @ L2[i, k, ])
        return L.mean(axis=1).prod().logabs

    def _Newton_bl_update(self, L1, L2, bl):
        '''Calculate the Newton's method update on the length of a branch.

        Parameters
        ----------
        L1, L2 : LogArray
            Conditional likelihood profiles at the two flanking nodes.
        bl : numpy.ndarray-like
            Initial branch length.
        '''
        L = LogArray(np.ones((self.nsite, self.ASRVk)))
        dL = LogArray(np.ones((self.nsite, self.ASRVk)))
        ddL = LogArray(np.ones((self.nsite, self.ASRVk)))
        for i in range(self.nsite):
            for k in range(self.ASRVk):
                PI_i_L1_ik = self.PI[i,] * L1[i, k, ]
                P_ik_L2_ik = self._calc_P(i, bl, k) @ L2[i, k, ]
                L[i, k] = PI_i_L1_ik @ P_ik_L2_ik
                dL[i, k] = PI_i_L1_ik @ (self.ASRVr[k] * self.Q[i,]) @ P_ik_L2_ik
                ddL[i, k] = PI_i_L1_ik @ ((self.ASRVr[k]**2) * (self.Q[i,] @ self.Q[i,])) @ P_ik_L2_ik
        dllk = (dL.sum(axis=1) / L.sum(axis=1)).sum().convert()
        ddllk = (((ddL.sum(axis=1) * L.sum(axis=1)) - (dL.sum(axis=1)**2)) / (L.sum(axis=1)**2)).sum().convert()
        return -dllk / ddllk

    def _optimize_bl_node(self, L1, L2, blInit, nodeLabel=None, verbose=False):
        '''Optimize the length of a branch.

        Parameters
        ----------
        L1, L2 : LogArray
            Conditional likelihood profiles at the two flanking nodes.
        blInit : numpy.ndarray-like
            Initial branch length.
        nodeLabel : int, optional
            Node label; only needs to be given if verbose is True.
        verbose : bool, optional
            Should the optimizationn process be printed? The default is False.
        '''
        # test for monotonically decreasing likelihood curve
        # if it is, branch length is set to the lower bound
        testBl = np.power(10., np.linspace(-4, 1, num=6))
        testllk = np.array([self._calc_llk_bl(L1, L2, bl) for bl in testBl])
        if np.all(np.diff(testllk) <= 0.):
            llk = self._calc_llk_bl(L1, L2, self._minbl)
            if verbose:
                print(f'   Node {nodeLabel}, monotonically decreasing likelihood curve.')
                print(f'      bl set to lower bound ({self._minbl:.6f}), llk = {llk:.6f}')
            return self._minbl, llk

        # first attempt Newton's method
        # if the initial update is out of bounds or decreases the likelihood,
        #   resort to bounded optimization
        llkInit = self.llk
        niter = 0
        if verbose:
            print(f'   Node {nodeLabel}, initial bl = {blInit:.6f}, initial llk = {llkInit:.6f}')

        while niter < 10:
            bl = blInit + self._Newton_bl_update(L1, L2, blInit)
            llk = self._calc_llk_bl(L1, L2, bl)
            # Anomaly
            if bl < self._minbl or bl > self._maxbl or llk < llkInit:
                if verbose:
                    print(f"      Round {niter + 1} Newton's method failed; " +
                          "trying bounded optimization.")
                break
            # Successful optimization
            if llk > llkInit and llk < llkInit + self._tol:
                if verbose: print(f"      Round {niter + 1} Newton's method convergence: " +
                                  f'bl = {bl:.6f}, llk = {llk:.6f}')
                return bl, llk
            # Otherwise, successful update
            if verbose:
                print(f"      Round {niter + 1} Newton's method update: bl = {bl:.6f}, llk = {llk:.6f}")
            blInit = bl
            llkInit = llk
            niter += 1

        # Bounded optimization
        def objective(bl):
            return -self._calc_llk_bl(L1, L2, bl)
        res = optimize.minimize_scalar(objective, bounds=(self._minbl, self._maxbl),
                                       method='bounded')
        if verbose:
            print(f'      bl = {res.x:.6f}, llk = {-res.fun:.6f}')
        return res.x, -res.fun

    def _optimize_bl_one_round(self):
        # Perform one round of branch length optimization.
        for node in self.tree.postorder_internal_node_iter():
            if node != self.tree.seed_node:
                if node.parent_node != self.tree.seed_node:
                    allk = self._calc_subtree_llk(node)
                    abl = node.edge_length
                elif self.tree.root_state == 2:
                    for sister in self.tree.seed_node.child_node_iter():
                        if sister != node:
                            allk = sister.llk
                    abl = 2. * node.edge_length
                else:
                    allk = np.zeros((self.nsite, self.ASRVk, 20))
                    for sister in self.tree.seed_node.child_node_iter():
                        if sister != node:
                            allk += self._prune_branch(sister.llk, sister.edge_length)
                    abl = node.edge_length
                for child in node.child_node_iter():
                    L1 = LogArray(child.llk, isLog=True)
                    L2 = self._prune_branch(allk, abl)
                    for sister in node.child_node_iter():
                        if sister != child:
                            L2 += self._prune_branch(sister.llk, sister.edge_length)
                    L2 = LogArray(L2, isLog=True)
                    bl, llk = self._optimize_bl_node(L1, L2, child.edge_length)
                    child.edge_length = bl
                    self.llk = llk
            else:
                if self.tree.root_state == 2:
                    C1, C2 = node.child_nodes()
                    L1 = LogArray(C1.llk, isLog=True)
                    L2 = LogArray(C2.llk, isLog=True)
                    bl, llk = self._optimize_bl_node(L1, L2, 2. * C1.edge_length)
                    C1.edge_length, C2.edge_length = (bl / 2., bl / 2.)
                    self.llk = llk
                else:
                    for child in node.child_node_iter():
                        L1 = LogArray(child.llk, isLog=True)
                        L2 = np.zeros((self.nsite, self.ASRVk, 20))
                        for sister in node.child_node_iter():
                            if sister != child:
                                L2 += self._prune_branch(sister.llk, sister.edge_length)
                        L2 = LogArray(L2, isLog=True)
                        bl, llk = self._optimize_bl_node(L1, L2, child.edge_length)
                        child.edge_length = bl
                        self.llk = llk
            node.llk = np.zeros((self.nsite, self.ASRVk, 20))
            for child in node.child_node_iter():
                node.llk += self._prune_branch(child.llk, child.edge_length)

    def _optimize_bl(self, tol=None):
        # Optimize the branch lengths
        if tol is None:
            tol = self._tol
        for niter in range(10):
            llkInit = self.llk
            self._optimize_bl_one_round()
            print(f'   Sequential branch length optimization round {niter + 1}: llk = {self.llk:.6f}')
            if self.llk < llkInit + tol:
                return None
        print('   Sequential branch length optimization ended without reaching convergence.')

    def _optimize_ASRV(self, tol=None):
        # Optimize the ASRV shape parameter.
        if self.ASRVk == 1:
            print('   ASRV model turned off.')
            return None

        def objective(a):
            self._set_ASRV(self.ASRVk, a if a.ndim == 0 else a[0])
            return -self._calc_llk()

        res = optimize.minimize(objective, self.ASRVa, method='L-BFGS-B',
                                bounds=[(self._minval, None)], tol=tol)
        if res.success:
            print(f'   ASRV shape parameter optimization: a = {self.ASRVa:.4f}, llk = {self.llk:.6f}.')
            return None
        else:
            res = optimize.minimize_scalar(objective, bounds=(self._minval, 10), method='bounded')
            if res.success:
                print(f'   ASRV shape parameter optimization: a = {self.ASRVa:.4f}, llk = {self.llk:.6f}.')
                return None
            else:
                print('   ASRV shape parameter optimization failed.')

    def _optimize_protein_model(self, tol=None):
        # Optimize the DNA model and the fitness function.
        def objective(param):
            if self.DNAModel == 'HKY85':
                kappa = param[0]
                pC, pG, pT = param[1:4]
                self._set_DNAModel(self.DNAModel, [1, kappa, 1, 1, kappa, 1], [1, pC, pG, pT])
            elif self.DNAModel == 'GTR':
                rGA, rGC, rTA, rTC, rTG = param[0:5]
                pC, pG, pT = param[5:8]
                self._set_DNAModel(self.DNAModel, [1, rGA, rGC, rTA, rTC, rTG], [1, pC, pG, pT])
            if self.fitnessFuncType == 'logistic':
                if self.DNAModel == 'JC69':
                    self._set_fitness_function(self.fitnessFuncType, param[0:3])
                elif self.DNAModel == 'HKY85':
                    self._set_fitness_function(self.fitnessFuncType, param[4:7])
                elif self.DNAModel == 'GTR':
                    self._set_fitness_function(self.fitnessFuncType, param[8:11])
            return -self.update_model()

        initParam = []
        bounds = []
        if self.DNAModel == 'HKY85':
            kappa = self.DNAR[2, 0]
            pC, pG, pT = self.DNAPI[1:4] / self.DNAPI[0]
            initParam += [kappa, pC, pG, pT]
            bounds += [(self._minval, None)] * 4
        elif self.DNAModel == 'GTR':
            rGA, rGC, rTA, rTC, rTG = self.DNAR[[2, 2, 3, 3, 3], [0, 1, 0, 1, 2]] / self.DNAR[1, 0]
            pC, pG, pT = self.DNAPI[1:4] / self.DNAPI[0]
            initParam += [rGA, rGC, rTA, rTC, rTG, pC, pG, pT]
            bounds += [(self._minval, None)] * 8
        if self.fitnessFuncType == 'logistic':
            initParam += self.fitnessFuncParam.tolist()
            bounds += [(self._mingr, self._maxgr)] + [(self._minval, None)] + [(None, None)]

        res = optimize.minimize(objective, initParam, method='L-BFGS-B', bounds=bounds, tol=tol)
        if res.success:
            print(f'   DNA model and fitness function joint optimization: llk = {self.llk:.6f}')
            print('   DNAR:\n', self.DNAR[np.tril_indices(4, -1)])
            print('   DNAPI:\n', self.DNAPI)
            print(f'   {self.fitnessFuncType} parameters:\n', self.fitnessFuncParam)
        else:
            print('   DNA model and fitness function joint optimization failed.')


# help command
if '-h' in sys.argv or '--help' in sys.argv:
    print('''
      Required arguments

      -a (or --alignment) : str
          Path to the phylip-format alignment file.

      -t (or --tree) : str
          Path to the newick-format tree file.

      -p (or --phi) : str
          Path to the csv-format DMS data file.
          Not required when useSiteHomogeneousModel = 1.

      -o (or --outdir) : str
          Path to an existing directory in which results should be written.


      Optional output-related arguments

      --ASR : 1 or 0
          Should ASR be performed? 1: Yes, 0: No (default). If yes, a tree with
          node labels showing reconstructed ancestral states is also generated.


      Optional model parameters

      -d (or --DNAModel) : str
          DNA model; one of 'GTR' (default), 'JC69', and 'HKY85'.

      -f (or --fitnessFuncType) : str
          Fitness function type; currently only 'logistic' (default).

      -k (or --ASRVk) : int
          Number of categories for the gamma-distributed among-site rate
          variation model. Default is 4. Set to 1 to disable the model.

      --DNAR : str of 6 floats separted by space
          Initial DNA exchangeability matrix. Enter the 6 lower triangular
          elements row-by-row (in the order A, C, G, T). Default is all 1.
          Example: '1 2 3 4 5 6'.

      --DNAPI : str of 4 floats separated by space
          Initial nucleotide frequency vector (A, C, G, T). Default is all 0.25.
          Example: '0.1 0.2 0.3 0.4'.

      --fitnessFuncParam : str of variable number of floats separated by space
          Initial fitness function parameters. Default is (1, 1, 0) for `logistic'.
          Example: '3 10 -0.5'.

      --ASRVa : float
          Initial shape parameter for the among-site rate variation model.
          Default is 1.

      --AAR : str
          Path to a file containing the 190 lower triangular elements of a
          site-homogeneous amino acid exchangeability matrix. File must contain
          190 lines of numbers entered row-by-row. Standard amino acid ordering
          (A, R, N, D, ...) is followed. This option is used only if
          --useSiteHomogeneousModel = 1.

      --AAPI : str
          Path to a file containing the site-homogeneous amino acid equilibrium
          frequency vector (A, R, N, D, ...). This option is used only if
          --useSiteHomogeneousModel = 1.


      Optional control arguments (1 = Yes, 0 = No)

      --optProteinModel : 1 (default) or 0
          Should the DMS model parameters be optimized?

      --optBl: 1 (default) or 0
          Should the branch lengths be optimized?

      --optASRV : 1 (default) or 0
          Should the among-site rate variation shape parameter be optimized?

      --useSiteHomogeneousModel : 1 or 0 (default)
          Should a standard site-homogeneous model be used? If yes, enter the
          model through --AAR and --AAPI.
          ''')
    exit()

# parsing command line inputs
short_options = 'a:t:p:o:d:f:k:'
long_options = ['alignment=', 'tree=', 'phi=', 'outdir=',
                'ASR=', 'DNAModel=', 'fitnessFuncType=', 'ASRVk=',
                'DNAR=', 'DNAPI=', 'fitnessFuncParam=', 'ASRVa=',
                'optProteinModel=', 'optBl=', 'optASRV=',
                'useSiteHomogeneousModel=', 'AAR=', 'AAPI=']
arg_val_list = getopt.getopt(sys.argv[1:], short_options, long_options)[0]

# default parameter values

DNAModel = 'GTR'
fitnessFuncType = 'logistic'
ASRVk = 4
DNAR = None
DNAPI = None
fitnessFuncParam = None
ASRVa = 1.

performASR = False
optProteinModel = True
optBl = True
optASRV = True
useSiteHomogeneousModel = False

# input parameter values

for arg, val in arg_val_list:
    # required arguments
    if arg in ('-a', '--alignment'):
        alignmentFile = val
    if arg in ('-t', '--tree'):
        treeFile = val
    if arg in ('-p', '--phi'):
        DMSFile = val
    if arg in ('-o', '--outdir'):
        outdir = val
    # optional output-related arguments
    if arg == '--ASR':
        if val == '1':
            performASR = True
    # optional model parameters
    if arg in ('-d', '--DNAModel'):
        DNAModel = val
    if arg in ('-f', '--fitnessFuncType'):
        fitnessFuncType = val
    if arg in ('-k', '--ASRVk'):
        ASRVk = int(val)
    if arg == '--DNAR':
        DNAR = [float(i) for i in val.split()]
    if arg == '--DNAPI':
        DNAPI = [float(i) for i in val.split()]
    if arg == '--fitnessFuncParam':
        fitnessFuncParam = [float(i) for i in val.split()]
    if arg == '--ASRVa':
        ASRVa = float(val)
    if arg == '--AAR':
        AAR = np.genfromtxt(val)
    if arg == '--AAPI':
        AAPI = np.genfromtxt(val)
    # optional control arguments
    if arg == '--optProteinModel':
        if val == '0':
            optProteinModel = False
    if arg == '--optBl':
        if val == '0':
            optBl = False
    if arg == '--optASRV':
        if val == '0':
            optASRV = False
    if arg == '--useSiteHomogeneousModel':
        if val == '1':
            useSiteHomogeneousModel = True
            optProteinModel = False
            DMSFile = ''

# initialization

model = DMSPhyloAA(alignmentFile, treeFile, DMSFile,
                   DNAModel, fitnessFuncType, ASRVk,
                   DNAR, DNAPI, fitnessFuncParam, ASRVa)

if useSiteHomogeneousModel:
    model._derive_protein_model(R=np.tile(AAR, (model.nsite, 1)),
                                PI=np.tile(AAPI, (model.nsite, 1)))
    model._calc_llk()

# optimization

if optProteinModel or optBl or optASRV:
    model.optimize(optProteinModel, optBl, optASRV)

# output

model.export(outdir)

if performASR:
    os.mkdir(outdir + 'ASR/')
    model.export_ASR(outdir + 'ASR/')
    model.export_state_tree(outdir)
