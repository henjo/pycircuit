# Based on m_matrix.h from Gnucap, Al Davis

# doctest not supported yet

from numpy import arange, array, zeros

class Bsmatrix:
  """Sparse matrix, bordered block with spikes
  
  declare a matrix:
  a=Bsmatrix(4)
  tell it what positions may be used
  a.iwant(0,1)
  a.iwant(1,2)
  a.iwant(0,3)
  setup the structure
  a.allocate()

  """
  def __init__(self,size=None,T=None):
    self.changed=[] #list with booleans
    self.data=None
    self.lownode=arange(size)
    self.row0=self.lownode-1  #index of position 0 of every row
    self.col0=self.lownode+1 #index of position 0 of every col
    self.diag=self.lownode # index of diagonal elements
    self.nzcount=0 # count of non-zero elements
    self.size=size
    self.min_pivot=0. # minimum pivot value
    
  def __getitem__(self, rc):
    row=rc[0]
    col=rc[1]
    assert(max(row,col)<=self.size)
    assert(min(row,col)>=0)
    if min(row,col)<self.lownode[max(row,col)]:
      return 0. #a "sparse 0"
    else:
      index=self.getindex(row,col)
      return self.data[index]

  def __setitem__(self, rc, val):
    row=rc[0]
    col=rc[1]
    assert(max(row,col)<=self.size)
    assert(min(row,col)>=0)
    assert(min(row,col)>=self.lownode[max(row,col)])
    index=self.getindex(row,col)
    self.data[index]=val
    return None

  def getindex(self,row,col):
    """Returns index of 1-dim data storage correponsing to 2-dim input

    The actual data is stored in a one dimensional array (self.data) 
    and getindex converts the 2-dimensional representation with row and 
    column, to the corresponding index in the 1-dim data array.

    >>> a=Bsmatrix(4)
    >>> a.iwant[0,1]
    >>> a.iwant[1,2]
    >>> a.iwant[0,3]
    >>> a.allocate()
    >>> a.getindex(2,1)
    4.0
    >>> a.getindex(2,2)
    5.0
    >>> a.getindex(1,2)
    6.0
    >>> a.getindex(3,0)
    7.0
    >>> a.getindex(0,2) #not part of sparse matrix
    7.0

    

    """
    if row>col: #below diagonal
      return self.row0[row]+col
    else: #on or above diagonal
      return self.col0[col]-row

  def is_changed(self,n):
    return self.changed[n]
  def set_changed(self,n, x = True):
    self.changed[n] = x
    return None

  def set_min_pivot(self,x):
    self.min_pivot = x
    return None
  
  def size(self):
    return self.size

  def new_subdot(self,i,j):
    """Solves for L/U elements by subtracting a dot-product

    Solves equation A(i.j)=L(i,:)*U(:,j) for L(i,j) or U(i,j) 
    """
    max_low_node=max(self.lownode[i],self.lownode[j])
    if i>j:# below diagonal
      min_node=j #first diagonal
      diagonal=self[j,j] #lower diagonal
    else:
      min_node=i
      diagonal=1
    for k in range(max_low_node,min_node-1):
      self[i,j] -= self[i,k]*self[k,j]
    self[i,j] /= diagonal
    return None

  def iwant(self, node1, node2):
    """Adapt self.lownode to matrix shape

    >>> a=Bsmatrix(4)
    >>> a.iwant[0,1]
    >>> a.iwant[1,2]
    >>> a.iwant[0,3]
    >>> a.allocate()
    >>> a.lownode
    array([0, 0, 1, 0])

    node 0 has a connection downto 0
    node 1 has a connection downto 0
    node 2 has a connection downto 1
    node 3 has a connection downto 0

    resulting structure (of a not so sparse matrix) : (used = X, unused=-)
    X X - X
    X X X X
    - X X X
    X X X X
    """
    assert(node1 < self.size)
    assert(node2 < self.size)
    if (node1 < 0 or node2 < 0):
      # negative is invalid, not used but still may be in a node list
      pass
    elif (node1 < self.lownode[node2]):
      self.lownode[node2]=node1
    elif (node2 < self.lownode[node1]):
      self.lownode[node1]=node2
      return None


  def allocate(self):
    """
    Setup based on matrix shape
    Now we update indices as well
    We calculate nzcount and create the data vector

    >>> a=Bsmatrix(4)
    >>> a.allocate()
    >>> a.data
    array([ 0.,  0.,  0.,  0.])
    >>> a.lownode=array([0, 0, 1, 2])
    >>> a.allocate()
    >>> a.data
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

    """
    low_diff=arange(self.size)-self.lownode
    self.nzcount=zeros(self.size)
    self.nzcount=2*(low_diff)+1
    self.nzcount=self.nzcount.cumsum()
    self.setzero() #makes sure all elements in space are zero
    self.diag=zeros(self.size)
    for i in xrange(self.size-1):
      self.diag[i+1]=self.diag[i]+low_diff[i]+low_diff[i+1]+1
    #self.diag=self.diag-1 #modification since matrix start@1
    self.row0=self.diag-arange(self.size)
    self.col0=self.diag+arange(self.size)
    return None

  def setzero(self):
    """
    setzero: makes sure all values are zero
    """
    self.data=zeros(self.nzcount[self.size-1])
    return None

  def dezero(self,offset):
    """
    This function adds and offset to the diagonal
    Although unlikely, the result might of course be zero
    To be sure we would have to check all diagonal elements
    for zero after the addition.

    >>> a=Bsmatrix(4)
    >>> a.allocate()
    >>> a.dezero(1)
    >>> a.dense()
    array([
    [ 1.,  0.,  0.,  0.],
    [ 0.,  1.,  0.,  0.],
    [ 0.,  0.,  1.,  0.],
    [ 0.,  0.,  0.,  1.]])

    """
    for i in xrange(self.size):
      self[i,i] += offset
    return None

  def lu_decomp(self):
    """LU-decomposition with normalised L-matrix
    Matrix self(=L*U) -> self=L+U, changed in place
    """
    for i in range(self.size):
      for j in range(self.lownode[i],i): #avoids diagonal
        self.new_subdot(i,j)
        self.new_subdot(j,i)
      self.new_subdot(i,i)
    return None

  def transpose(self):
    """Transpose matrix in place
    
    Swap all *valid* elements (i,j)->(j,i) not on diagonal.
    """
    for i in range(self.size):
      for j in range(self.lownode[i],i): #avoids diagonal
        temp=self[i,j]
        self[i,j]=self[j,i]
        self[j,i]=temp
      return None


  def fbsub(self,vec):
    """
    Forward backward substitution
    vec is right side vector and will be substituted in place
    """
    for i in xrange(self.size):
      for k in xrange(i):
        vec[i] -= vec[k]*self[i,k]
    for i in xrange(self.size):
      for k in xrange(self.size-1,i,-1):
        vec[i] -= vec[k]*self[i,k]
        vec[i] /= self[i,i]
    return None

  def dense(self):
    """
    Returns a dense array version of the sparse matrix
    """
    dense_arr=zeros((self.size,self.size))
    for i in xrange(self.size):
      for j in xrange(self.lownode[i],i):
        dense_arr[i,j]=self[i,j]
        dense_arr[j,i]=self[j,i]
      dense_arr[i,i]=self[i,i]
    return dense_arr

        
##Note! lu_decomp, transpose and dense are very similar
##in traversing the valid values of the sparse matrix
##Could be a good idea to have a traversal function, like for element in
##Is that called an iterator?

# >>> a=Bsmatrix(4)
# >>> a.iwant[0,1]
# >>> a.iwant[1,2]
# >>> a.iwant[0,3]
# >>> a.allocate()
# >>> a.lownode

    
