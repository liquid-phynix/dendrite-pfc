def read_array(fname):
    ret = None
    with open(fname, mode = 'rb') as fhandle:
        meta = fhandle.readline()
        if meta[0] != 35: raise ValueError('wrong metadata formats')
        meta = dict(filter(lambda x: len(x) == 2, [kv.split('=') for kv in meta.decode().strip('# ').split()]))
        ret = fromfile(fhandle, dtype = meta['dtype']).reshape((int(meta['rows']), int(meta['cols'])))
    return ret
def save_array(fname, arr):
    import sys
    order = sys.byteorder
    dtype = arr.dtype.newbyteorder('<' if order == 'little' else '>')
    with open(fname, mode = 'wb') as fhandle:
        fhandle.write(('# rows=%d cols=%d dtype=%s\n' % (arr.shape[0], arr.shape[1], dtype.str)).encode())
        arr.tofile(fhandle)

def plot_cont(num):
    points = read_array('points_' + num + '.txt')
    edges = read_array('endpoint_indices_' + num + '.txt')
    comu = points[edges]
    plot([comu[:,0,0],comu[:,1,0]], [comu[:,0,1],comu[:,1,1]], 'black')
    #    plot([comu[:,0,0]+0.5,comu[:,1,0]+0.5], [comu[:,0,1]+0.5,comu[:,1,1]+0.5], 'black')
