#!/usr/bin/env python

import argparse as ap
import tables
import numpy as np
from scipy.sparse import csr_matrix
import sys
import logging

def toString(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s

    if isinstance(s, bytes):  # or isinstance(s, np.bytes_):
        if sys.version_info[0] == 2:
            return str(s)
        return s.decode('ascii')

    if isinstance(s, list):
        return [toString(x) for x in s]

    if isinstance(s, np.ndarray):
        return s.astype(str)

    return s


def loadH5subset(filename, includechroms=None):
    '''
    loadH5(filename, includechroms=None, csr=True)

    loads an *.h5 hic matrix as created by hicexplorer

    :param filename:        name of the *.h5 file containing the matrix
    :param includechroms:   list of chromosomes to include in the returned objects
                            if not given all chromosomes in the *.h5 file are included
    :param csr:             if True returns a csr_matrix object else a full numpy.array

    :return:                csr_matrix containing the data in the matrix
    '''
    with tables.open_file(filename) as f:
        parts = {}
        try:
            for matrix_part in ('data', 'indices', 'indptr', 'shape'):
                parts[matrix_part] = getattr(f.root.matrix, matrix_part).read()
        except Exception:
            logging.info('No h5 file. Please check parameters concerning the file type!')
            exit(1)

        matrix = csr_matrix(tuple([parts['data'], parts['indices'], parts['indptr']]),
                            shape=parts['shape'], dtype = int)

        intervals = {}
        for interval_part in ('chr_list', 'start_list', 'end_list', 'extra_list'):
            if toString(interval_part) == toString('chr_list'):
                chrom_list = getattr(f.root.intervals, interval_part).read()
                intervals[interval_part] = toString(chrom_list)
            else:
                intervals[interval_part] = getattr(f.root.intervals, interval_part).read()

        cut_intervals = list(
            zip(intervals['chr_list'], intervals['start_list'], intervals['end_list'], intervals['extra_list']))

        assert len(cut_intervals) == matrix.shape[0], \
            "Error loading matrix. Length of bin intervals ({}) is different than the " \
            "size of the matrix ({})".format(len(cut_intervals), matrix.shape[0])

        # compute index array and chromosome list
        inds, chr_list, chroms = [], [], set()
        for i, (chr, start, end, extra) in enumerate(cut_intervals):
            if chr not in chroms:
                chroms.add(chr)
                inds.append(i)
                chr_list.append(chr)

        # if includechroms is given we filter the output for the chromosomes listed
        # and recompute indices of chromosome boundaries in the resulting matrix
        if includechroms:
            includechroms = set(includechroms)

            intervals = {k: [] for k in ['chr_list', 'start_list', 'end_list', 'extra_list']}
            for i, vals in enumerate(cut_intervals):
                if vals[0] in includechroms:
                    for k, v in zip(['chr_list', 'start_list', 'end_list', 'extra_list'], vals):
                        intervals[k].append(v)

            filterinds, filterchrs = [], []
            for i, chr in zip(range(len(inds)), chr_list):
                if chr in includechroms:
                    filterinds.append([inds[i], inds[i + 1] if i + 1 != len(inds) else matrix.shape[0]])
                    filterchrs.append(chr)

            matrixinds = np.zeros(shape=matrix.shape[0], dtype=bool)
            ncuts, tmpe = [], 0
            for s, e in filterinds:
                matrixinds[s: e] = True

                if s == tmpe:
                    ncuts.append(s)
                    tmpe = e

                else:
                    ncuts.append(tmpe)
                    tmpe = e - s + tmpe

            matrix = matrix[matrixinds, :][:, matrixinds]

    return matrix, intervals

def writeH5subset(matrix, intervals, outfile):
    '''
    writes a given matrix to *.h5 format as specified by the HiCExplorer tools

    :param matrix:      matrix in either np.ndarray or csr_matrix format
    :param inds:        indices denoting the start of each chromosome in the contact matrix
    :param chrlist:     list of chromosomes included in the contact matrix
    :param outfile:     file to write the matrix to

    :return:            None
    '''
    # initialize output file
    filter = tables.Filters(complevel = 5, complib = 'blosc')
    with tables.open_file(outfile, mode = 'w', title = 'Normalized Matrix') as h5file:
        matrix_group = h5file.create_group('/', 'matrix', )

        # save csr components
        for attr in ['data', 'indices', 'indptr', 'shape']:
            arr = np.array(getattr(matrix, attr))
            atom = tables.Atom.from_dtype(arr.dtype)
            ds = h5file.create_carray(matrix_group, attr, atom, shape = arr.shape, filters = filter)
            ds[:] = arr

        interval_group = h5file.create_group('/', 'intervals', )
        # save additional information]
        for k, v in intervals.items():
            arr = np.array(v)
            atom = tables.Atom.from_dtype(arr.dtype)
            ds = h5file.create_carray(interval_group, k, atom, shape = arr.shape, filters = filter)
            ds[:] = arr

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-m', '--matrix', required = True,
                    help = '*.h5 file to be processed')
parser.add_argument('--includeChroms', nargs = '+',
                    help = 'space-separated list of chromosomes that should be extracted from the given file')
parser.add_argument('-o', '--outfile', required = True,
                    help = 'file to which the resulting subset matrix should be written')
args = parser.parse_args()

writeH5subset(*loadH5subset(args.matrix, includechroms = args.includeChroms), args.outfile)
