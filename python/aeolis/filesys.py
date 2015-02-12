import os
import struct
import numpy as np
import pandas as pd

def load_dimensions(fpath):
    '''Load model dimensions from dims.out file

    Parameters
    ----------
    fpath : string
        Path to model location

    Returns
    -------
    dict
        Model dimensions, space and time axes
    '''

    if not os.path.isdir(fpath):
        fpath = os.path.split(fpath)[0]
        
    nx, dx, nt, dt, nf, nl, nt_out, dt_out = \
        read_fortran(os.path.join(fpath, 'dims.out'))

    # integers
    nx = int(nx)+1
    nt = int(nt)
    nf = int(nf)
    nl = int(nl)
    nt_out = int(nt_out)
    
    # axes
    ax_t = np.arange(nt) * dt
    ax_t_out = np.arange(nt_out) * dt_out
    ax_x = np.arange(nx) * dx
    
    return {'nx':nx,
            'dx':dx,
            'nt':nt,
            'dt':dt,
            'nf':nf,
            'nl':nl,
            'nt_out':nt_out,
            'dt_out':dt_out,
            'ax_t':ax_t,
            'ax_t_out':ax_t_out,
            'ax_x':ax_x}


def get_dims(fname):
    '''Determine file dimensions

    Parameters
    ----------
    fname : string
        Path to output file

    Returns
    -------
    tuple
        Tuple with file dimensions
    '''

    fpath, fname = os.path.split(fname)
    fname = fname.split('.')[0]
    d = load_dimensions(fpath)

    if fname in ['mass']:
        return (-1, d['nx'], d['nl'], d['nf'])
    elif fname in ['d10', 'd50', 'd90']:
        return (-1, d['nx'], d['nl'])
    elif fname in ['Cu', 'Ct', 'uth', 'supply']:
        return (-1, d['nx'], d['nf'])
    elif fname in ['x', 'z', 'moist_map']:
        return (-1, d['nx'])
    elif fname in ['rho', 'dist']:
        return (-1, d['nf'])
    else:
        return [-1]

    
def guess_dims(dims, data, start=0, stop=np.inf, step=1):
    '''Guess dimensions of data read from output file

    Parameters
    ----------
    dims : dict
        Dimensions read from dims.out file
    data : np.ndarray
        Data read from output file
    start : int, optional
        Start position used for data read
    stop : int, optional
        Stop position used for data read
    step : int, optional
        Step or stride used for data read

    Returns
    -------
    tuple
        Tuple with guessed dimensions of read data
    '''
    
    if np.isfinite(stop):
        n = int((stop - start) / step) + 1
    else:
        n = int((dims['nt_out'] - start) / step) + 1
        
    if check_dim(n, np.prod(data.shape) / dims['nx'] / dims['nf'] / dims['nl']) and \
       dims['nl'] * dims['nf'] > 1:
        return (-1, dims['nx'], dims['nl'], dims['nf'])
    elif check_dim(n, np.prod(data.shape) / dims['nx'] / dims['nl']) and dims['nl'] > 1:
        return (-1, dims['nx'], dims['nl'])
    elif check_dim(n, np.prod(data.shape) / dims['nx'] / dims['nf']) and dims['nf'] > 1:
        return (-1, dims['nx'], dims['nf'])
    elif check_dim(n, np.prod(data.shape) / dims['nx']):
        return (-1, dims['nx'])
    else:
        print np.prod(data.shape), dims['nx'], dims['nf'], dims['nl'], n
        return [-1]

    
def load_dataframe(fname, shape=None, start=0, stop=np.inf, step=1, verbose=False):
    '''Load data from output file in Pandas dataframe or equivalent

    Parameters
    ----------
    fname : string
        Path to output file
    shape : tuple or list, optional
        Shape of data in file
    start : int, optional
        Start block position for read
    stop : int, optional
        Stop block position for read
    step : int, optional
        Step or stride blocks for read
    verbose : bool, optional
        Flag to enable process output

    Returns
    -------
    pandas.DataFrame, pandas.Panel or pandas.Panel4D
        Pandas object containing data from output file

    Examples
    --------
    >>> load_dataframe('mass.out', stop=10) # read first 10 blocks
    >>> load_dataframe('supply.out', step=10) # read every 10th block
    >>> load_dataframe('supply.out', start=100, step=2, stop=200) # read every even block from 100th to 200th
    '''

    if shape is None:
        shape = get_dims(fname)

    dims = load_dimensions(fname)
    data = read_fortran(fname, shape=shape,
                        start=start, stop=stop, step=step, verbose=verbose)
        
    #ix = pd.TimedeltaIndex(start=0, periods=z.shape[0]-1, freq='H')
    ix = pd.DatetimeIndex(start=0, periods=data.shape[0], freq='%dS' % round(dims['dt_out'] * step))
    
    if len(shape) == 1:
        return pd.DataFrame(data, index=ix)
    elif len(shape) == 2:
        return pd.DataFrame(data, index=ix, columns=dims['ax_x'])
    elif len(shape) == 3:
        return pd.Panel(data, items=ix, major_axis=dims['ax_x'])
    elif len(shape) == 4:
        return pd.Panel4D(data, labels=ix, items=dims['ax_x'])
    else:
        raise ValueError('Unsupported dimension count [%d]' % len(shape))
    
    
def read_fortran(fname, shape=None, start=0, stop=np.inf, step=1, verbose=False):
    '''Load data from output file in numpy array

    Parameters
    ----------
    fname : string
        Path to output file
    shape : tuple or list, optional
        Shape of data in file
    start : int, optional
        Start block position for read
    stop : int, optional
        Stop block position for read
    step : int, optional
        Step or stride blocks for read
    verbose : bool, optional
        Flag to enable process output

    Returns
    -------
    np.ndarray
        numpy array containing data from output file
    '''
    
    s = 0
    A = []
    if os.path.exists(fname):
        with open(fname, 'rb') as fp:
            while True:
                try:
                    n = struct.unpack('I', fp.read(4))[0]
                    s += 1
                    
                    if (len(A) == 0 and s <= start) or (len(A) > 0 and s < step):
                        fp.seek(n + 4, 1)
                        continue
                    
                    if verbose:
                        print 'Read block #%d of size %d' % (len(A)+1, n)
                            
                    contents = fp.read(n)
                    if n == 4:
                        Ai = struct.unpack('i', contents)
                    else:
                        Ai = struct.unpack('%dd' % (n/8), contents)
                        
                    fp.seek(4, 1) # redundant
                        
                    if len(Ai) == 1:
                        A.append(Ai[0])
                    else:
                        A.append(np.asarray(Ai))
                            
                    if len(A) >= stop / step:
                        if verbose:
                            print 'Read maximum number of blocks. Abort.'
                                
                        break
                                        
                    s = 0
                                        
                except:
                    break
                    
    A = np.asarray(A)
                    
    if shape is not None:
        A = A.reshape(shape)
            
    return A
                    
                    
def check_dim(n, x):
    '''Check if guessed dimensions are valid'''
    
    if not is_int(x):
        return False
    elif abs(x - n) > 10:
        return False
    else:
        return True
    
        
def is_int(x):
    '''Check if variable is an integer'''
    
    return x - round(x) == 0.
    
