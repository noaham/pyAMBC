def aorn(a,n):
    if a == 0:
        return n
    else:
        return a

def is_partition(L):
    # checks if a list of integers is a partition (i.e. is it weakly decreasing)
    for i in range(len(L))[:-1]:
        if L[i+1] > L[i]:
            return False
    return True

def is_subpartition(mu,lam):
    # checks if the diagram of mu fits in the diagram of lambda
    if (is_partition(lam) and is_partition(mu)):
        if len(mu) > len(lam):
            return False
        else:
            for i in range(len(mu)):
                if mu[i] > lam[i]:
                    return False
        return True
    else:
        return False

def strip_zeros(mu):
    # strips all zeros from a list
    return [x for x in mu if x != 0]

def strip_empties(T):
    # strips trailing []'s from end of list.
    if T == []:
        return T
    elif T[-1] == []:
        return strip_empties(T[:-1])
    return T
        
    
def transpose_partition(lam):
    # finds the transpose of a partition
    if lam == []:
        return []
    return [len([j for j in lam if j > i]) for i in range(lam[0])]

def outer_addrem(lam):
    # produces a list of outer addable and removable boxes of the partition (i.e. boxes
    # that can be added or removed on the southeast edge)
    addables = []
    removables = []
    for i in range(1,len(lam)+1):
        j = lam[i-1]+1
        newpart = lam[:]
        newpart[i-1] += 1
        if is_partition(newpart):
            addables.append((i,j))
        j = lam[i-1]
        newpart = lam[:]
        newpart[i-1] -= 1
        if is_partition(newpart):
            removables.append((i,j))
    addables += [(len(lam)+1,1)]
    return [addables,removables]
    
class skewshape(object):
    # defines the class of skew shapes
    def __init__(self,lam, mu=[]):
        self.outer_partition = strip_zeros(lam)
        self.inner_partition = strip_zeros(mu)
        # the below is mu padded with zero rows to make it the same length as lamda
        self.inner_partition_ext = mu + [0]*(len(lam)-len(mu))
        if not is_subpartition(mu,lam):
            raise ValueError("Lists must define partitions that give a valid skewshape")
        self.row_sizes = [self.outer_partition[i] - self.inner_partition_ext[i] for i in range(len(self.outer_partition))]
        self.numboxes = sum(self.row_sizes) # sum(self.outer_partition) - sum(self.inner_partition)
        self.boxes = []
        for i in range(1,len(self.outer_partition)+1):
            self.boxes += [(i,j) for j in range(self.inner_partition_ext[i-1]+1,self.outer_partition[i-1]+1)]
        self.outer_addables = outer_addrem(self.outer_partition)[0]
        self.outer_removables = [box for box in outer_addrem(self.outer_partition)[1] if box in self.boxes]
        self.inner_addables = outer_addrem(self.inner_partition)[1]
        self.inner_removables = [box for box in outer_addrem(self.inner_partition)[0] if box in self.boxes]

    def __repr__(self):
        return str(self.outer_partition) + ' / ' + str(self.inner_partition)

    def __eq__(self,other):
        return (self.inner_partition == other.inner_partition) and (self.outer_partition == other.outer_partition)
        
    def diagram(self):
        # prints out a diagram representing the skew shape
        repstr = ''
        lam = self.outer_partition[:]
        mu = self.inner_partition[:]
        mu += [0]*(len(lam)-len(mu))
        for i in range(len(self.outer_partition)):
            repstr += '. '*mu[i] + 'x '*(lam[i]-mu[i]) + '\n'
        return repstr

    def transpose(self):
        # transposes a skew shape
        inner = transpose_partition(self.inner_partition)
        outer = transpose_partition(self.outer_partition)
        return skewshape(outer,inner)

    def is_box(self,i,j):
        # checks if position (i,j) is a box contained in the diagram of the skew shape
        if i < 1 or j < 1:
            return False
        elif len(self.outer_partition) < i:
            return False
        elif self.outer_partition[i-1] < j:
            return False
        elif self.inner_partition_ext[i-1] >= j:
            return False
        return True

    def __contains__(self,box):
        return self.is_box(box[0],box[1])

    def rembox(self,i,j, type=None):
        # remove the (i,j) box
        if (i,j) in self.inner_removables and type != 'outer':
            newmu = self.inner_partition[:]
            if len(newmu) < i:
                newmu = newmu + [1]
            else:
                newmu[i-1] += 1
            return skewshape(self.outer_partition,newmu)
        elif (i,j) in self.outer_removables and type != 'inner':
            newlam = self.outer_partition[:]
            newlam[i-1] -= 1
            return skewshape(newlam,self.inner_partition)
        else:
            raise ValueError("must be a removable box (of appropriate type)")

    def addbox(self,i,j):
        # adds a box (i,j)
        if (i,j) in self.outer_addables:
            newouter = self.outer_partition[:]
            if len(newouter) < i:
                newouter.append(1)
            else:
                newouter[i-1] += 1
            return skewshape(newouter,self.inner_partition)
        elif (i,j) in self.inner_addables:
            newmu = self.inner_partition[:]
            newmu[i-1] -= 1
            return skewshape(self.outer_partition,newmu)
        else:
            raise ValueError("must be an addable box")

    def __add__(self,other):
        # defines the sum (or union). Here other must be a shape that extends self.
        if self.outer_partition == other.inner_partition:
            mu = self.inner_partition[:]
            lam = other.outer_partition[:] 
            return skewshape(lam,mu)
        else:
            raise ValueError("shapes must extend each other")



class tableau(object):
    # class of a tableau in a skew shape. sh is a skewshape and T is a list of lists,
    # representing the rows (including empty rows).
    def __init__(self,sh,T):
        self.shape = sh
        self.T = T

        self.rows = [['.']*(self.shape.inner_partition_ext[i])+self.T[i] for i in range(len(self.T))]
        self.entries = [a for row in self.T for a in row]
        if not self.has_hole() and self.entries != []:
            self.max_entry = max(self.entries)

        if self.shape.row_sizes != [len(row) for row in T]:
            raise ValueError("T must define an entry for every available box in each row.")
        if self.is_semisemistd():
            self.boxes_ordered = self.order_boxes()

    def lt(self,a,b):
        # a and b are two boxes (in the form (i,j)). This returns true if a is less than
        # b, in the sense that either the entry in a is < the entry in b, or a occures to
        # the left of b.
        if self.is_semisemistd():
            if self.entry(a[0],a[1]) == self.entry(b[0],b[1]):
                if a[1] == b[1]:
                    return a[0] < b[0]
                return a[1] < b[1]
            else:
                return self.entry(a[0],a[1]) < self.entry(b[0],b[1])
        else:
            raise ValueError("tableau must be semistandard for lt to work")

    def order_boxes(self):
        boxes = self.shape.boxes[:]
        sorted = False
        while not sorted:
            sorted = True
            for i in range(len(boxes)-1):
                if not self.lt(boxes[i],boxes[i+1]):
                    a = boxes[i]
                    b = boxes[i+1]
                    boxes[i] = b
                    boxes[i+1] = a
                    sorted = False
        return boxes
    
    def col(self,j):
        # returns the jth column
        return [row[j-1] for row in self.rows if len(row) >= j]

    def entry(self,i,j):
        if self.shape.is_box(i,j):
            return self.rows[i-1][j-1]
        else:
            raise ValueError("must be a valid box")

    def __getitem__(self,box):
        if box in self.shape:
            return self.rows[box[0]-1][box[1]-1]
        else:
            raise ValueError("must be a valid box")
    
    def __eq__(self,other):
        return (self.T == other.T) and (self.shape == other.shape)

    def __repr__(self):
        if self.rows == []:
            return "empty"
        rep_str = ""
        colwidths = ['.' for i in range(len(self.rows[0]))]
        for j in range(len(self.rows[0])):
            max_digits = 1
            for entry in self.col(j+1):
                digits = len(str(entry))
                max_digits = max([max_digits,digits])
            colwidths[j] = max_digits
        for row in self.rows:
            row_str = ""
            for i in range(len(row)):
                space_buffer = colwidths[i] - len(str(row[i]))
                row_str += str(row[i])+" "*(space_buffer+1)
            rep_str += row_str +"\n"
        return rep_str[:-2]

    def __mul__(self,other):
        # multiplication in the plactic monoid
        if self.shape.outer_partition == []:
            return other
        if other.shape.outer_partition == []:
            return self
        shift = self.shape.outer_partition[0]
        newsh_outer = [a + shift for a in other.shape.outer_partition] + self.shape.outer_partition
        newsh_inner = [a + shift for a in other.shape.inner_partition_ext] + self.shape.inner_partition
        newsh = skewshape(newsh_outer,newsh_inner)
        return tableau(newsh, other.T + self.T).rect()
        

    def has_hole(self):
        # determines whether or not T has a hole. Box is given by self.hole
        i=0
        for row in self.rows:
            i+=1
            if 'x' in row:
                self.hole = (i,row.index('x')+1)
                return True
        return False
    
    def transpose(self):
        shape_trans = self.shape.transpose()
        T_trans = [self.col(i)[shape_trans.inner_partition_ext[i-1]:] for i in range(1,self.outer_partition[0]+1)]
        return tableau(shape_trans,T_trans)

    def is_semistd(self):
        # checks if the tableau has weakly increasing rows and increasing columns
        if len(self.entries) < 2:
            return True
        if self.has_hole():
            return False
        for row in self.rows:
            for j in range(1,len(row)):
                if row[j-1] != '.':
                    if row[j-1] > row[j]:
                        return False
        for colu in [self.col(j) for j in range(1,len(self.rows[0])+1)]:
            for i in range(1,len(colu)):
                if colu[i-1] != '.':
                    if colu[i-1] >= colu[i]:
                        return False
        return True

    def is_semisemistd(self):
        # checks if the tableau has weakly increasing rows and weakly increasing columns
        if len(self.entries) < 2:
            return True
        if self.has_hole():
            return False
        for row in self.rows:
            for j in range(1,len(row)):
                if row[j-1] != '.':
                    if row[j-1] > row[j]:
                        return False
        for colu in [self.col(j) for j in range(1,len(self.rows[0])+1)]:
            for i in range(1,len(colu)):
                if colu[i-1] != '.':
                    if colu[i-1] > colu[i]:
                        return False
        return True

    def is_std(self):
        # checks if the tableau is standard
        return self.is_semistd() and (len(set(self.entries)) == len(self.entries))
    
    def rembox(self,i,j,type=None):
        # remove the (i,j) box and either: create a hole if it is not a removable
        # node, or delete the box and change the shape if it isn't
        if self.shape.is_box(i,j):
            newT = [a[:] for a in self.T]
            if (i,j) in self.shape.inner_removables and type != 'outer':
                newT[i-1] = newT[i-1][1:]
                newsh = self.shape.rembox(i,j,type)
            elif (i,j) in self.shape.outer_removables and type != 'inner':
                newT[i-1] = newT[i-1][:-1]
                newsh = self.shape.rembox(i,j,type)
            else:
                newT[i-1][j-1-self.shape.inner_partition_ext[i-1]] = 'x'
                newsh = self.shape
            if (len(newT) == len(newsh.outer_partition) + 1) and newT[-1] == []:
                newT = newT[:-1]
            return tableau(newsh,newT)
        else:
            raise ValueError("must be a valid box")

    def addbox(self,i,j,x):
        # adds a box (i,j) filled with x
        if (i,j) in self.shape.outer_addables:
            newT = [a[:] for a in self.T]
            if len(newT) < i:
                newT = newT + [[x]]
            else:
                newT[i-1] = newT[i-1] + [x]
            return tableau(self.shape.addbox(i,j),newT)
        elif (i,j) in self.shape.inner_addables:
            newT = [a[:] for a in self.T]
            newT[i-1] = [x] + newT[i-1]
            return tableau(self.shape.addbox(i,j),newT)
        else:
            raise ValueError("must be an addable box")

    def __add__(self,other):
        # defines the sum (or union). Here other must be a shape that extends self.
        add_shape = self.shape + other.shape
        max_rows = max([len(self.T),len(other.T)])
        selfT_ext = self.T[:] +[[]]*(max_rows - len(self.T))
        othT_ext = other.T[:] +[[]]*(max_rows - len(other.T))
        newT = [selfT_ext[p] + othT_ext[p] for p in range(max_rows)]
        return tableau(add_shape,newT)
        
    def repl_entry(self,i,j,x):
        # replaces the entry in T in box (i,j) with x
        newT = [a[:] for a in self.T]
        newT[i-1][j-1-self.shape.inner_partition_ext[i-1]] = x
        return tableau(self.shape, newT)

    def restrict(self,p,q):
        # restricts a tableaux to entires in the range [p,q]. Must have weakly
        # increasing rows and columns.
        if self.is_semisemistd():
            Tres = self
            for box in [box for box in self.boxes_ordered[::-1] if self[box] > q]:
                Tres = Tres.rembox(box[0],box[1],'outer')
            for box in [box for box in self.boxes_ordered if self[box] < p]:
                Tres = Tres.rembox(box[0],box[1],'inner')
            return Tres
        else:
            raise ValueError("must have weakly increasing rows and columns")
        
    def rev_jdt(self):
        # returns the result of performing a single reverse jue de taquin slide
        # into the box into the hole marked by 'x'
        if self.has_hole():
            i = self.hole[0]
            j = self.hole[1]
            if self.hole in self.shape.outer_removables:
                return self.rembox(i,j,'outer')
            poss_slides = [box for box in [(i,j+1),(i+1,j)] if self.shape.is_box(box[0],box[1])]
            entries = [self.entry(box[0],box[1]) for box in poss_slides]
            ent = min(entries)
            if len(poss_slides) == 2 and entries[0] == entries[1]:
                sliding_box = (i+1,j)
            else:
                sliding_box = poss_slides[entries.index(min(entries))]
            newT = self.repl_entry(i,j,ent).repl_entry(sliding_box[0],sliding_box[1],'x')
            return newT
        else:
            raise ValueError("must have a hole to slide into")

    # def rev_slide(self,i,j,return_box=False):
    #     # result of doing full reverse slide into the inner addable (i,j)
    #     if (i,j) in self.shape.inner_addables:
    #         Tx = self.addbox(i,j,'x')
    #         while Tx.has_hole():
    #             Tx = Tx.rev_jdt()
    #         if return_box:
    #             big = self.shape.outer_partition[:]
    #             small = Tx.shape.outer_partition[:]
    #             box_shape = skewshape(big,small)
    #             return [Tx,box_shape]
    #         else:
    #             return Tx
    #     else:
    #         raise ValueError("must be an inner addable to slide into")

    def rev_slide(self,S,return_evac=False):
        # returns the result of reverse sliding self into the boxes determined by a
        # semistandard tableau S. The shape of self must extend the shape of S.
        if self.shape.inner_partition == S.shape.outer_partition:
            evac_boxes = [self.shape.outer_partition]
            boxes = reversed(S.boxes_ordered[:])
            Tx = self
            for box in boxes:
                i = box[0]
                j = box[1]
                Tx = Tx.addbox(i,j,'x')
                while Tx.has_hole():
                    Tx = Tx.rev_jdt()
                #box_shape = skewshape(big,small)
                evac_boxes = [Tx.shape.outer_partition] + evac_boxes
            if return_evac:
                #print(evac_boxes)
                return [Tx,chain_to_tableau(evac_boxes)]
            else:
                return Tx
        else:
            raise ValueError("shape of S must be extended by self")
            

    # def rect(self,return_evac=False):
    #     # if self.shape.inner_addables == []:
    #     #     return self
    #     newT = self
    #     num_addables = len(self.shape.inner_addables)
    #     outer = self.shape.outer_partition
    #     evac_chain = []
    #     while num_addables > 0:
    #         box = newT.shape.inner_addables[0]
    #         i = box[0]
    #         j = box[1]
    #         [newT,box] = newT.rev_slide(i,j,True)
    #         evac_chain = [box] + evac_chain
    #         num_addables = len(newT.shape.inner_addables)
    #     if return_evac:
    #         return [newT,chain_to_tableau(evac_chain)]
    #     return newT

    def rect(self,return_evac=False):
        Y = YamTab(skewshape(self.shape.inner_partition))
        return self.rev_slide(Y,return_evac)

    def fwd_jdt(self):
        # returns the result of performing a single forward jue de taquin slide
        # into the box into the hole marked by 'x'
        if self.has_hole():
            (i,j) = self.hole
            if self.hole in self.shape.inner_removables:
                return self.rembox(i,j,'inner')
            poss_slides = [box for box in [(i,j-1),(i-1,j)] if self.shape.is_box(box[0],box[1])]
            entries = [self.entry(box[0],box[1]) for box in poss_slides]
            ent = max(entries)
            if len(poss_slides) == 2 and entries[0] == entries[1]:
                sliding_box = (i-1,j)
            else:
                sliding_box = poss_slides[entries.index(max(entries))]
            newT = self.repl_entry(i,j,ent).repl_entry(sliding_box[0],sliding_box[1],'x')
            return newT
        else:
            raise ValueError("must have a hole to slide into")

    # def fwd_slide(self,i,j):
    #     # result of doing full forward slide into the outer addable (i,j)
    #     if (i,j) in self.shape.outer_addables:
    #         Tx = self.addbox(i,j,'x')
    #         while Tx.has_hole():
    #             Tx = Tx.fwd_jdt()
    #         return Tx
    #     else:
    #         raise ValueError("must be an outer addable to slide into")

    def fwd_slide(self,S,return_evac=False):
        # returns the result of forward sliding self into the boxes determined by a
        # semistandard tableau S. The shape of S must extend the shape of self.
        if self.shape.outer_partition == S.shape.inner_partition:
            evac_boxes = [self.shape.inner_partition]
            boxes = S.boxes_ordered[:]
            Tx = self
            for box in boxes:
                i = box[0]
                j = box[1]
                Tx = Tx.addbox(i,j,'x')
                while Tx.has_hole():
                    Tx = Tx.fwd_jdt()
                evac_boxes = evac_boxes + [Tx.shape.inner_partition]
            if return_evac:
                return [Tx,chain_to_tableau(evac_boxes)]
            else:
                return Tx
        else:
            raise ValueError("shape of S must extend the shape of self")
    
    def schutzenberger(self,n=None):
        # applies the schutzenberger involution
        if self.shape.boxes == []:
            return self
        if n is None:
            n = self.max_entry
        # if len(self.T)==1 and len(self.T[0])==1:
        #     return tableau([[n+1-self.T[0][0]]])
        if self.is_semistd():
            # first find smallest entry
            box = self.boxes_ordered[0]
            new_outer = self.shape.inner_partition[:]
            if box[0] > len(new_outer):
                new_outer += [1]
            else:
                new_outer[box[0]-1] += 1
            #box_sh = skewshape(new_outer,self.shape.inner_partition)
            m = self.entry(box[0],box[1])
            # now delete this box, make one slide and find the outer box
            # that was emptied
            box_tab = chain_to_tableau([self.shape.inner_partition,new_outer])
            [one_step,box_evac] = self.rembox(box[0],box[1],'inner').rev_slide(box_tab,True)
            box_evac = box_evac.repl_entry(box_evac.shape.boxes[0][0],box_evac.shape.boxes[0][1],n+1-m)
            # big = self.shape.outer_partition[:]
            # small = one_step.shape.outer_partition[:]
            # if len(small) < len(big):
            #     small = small + [0]
            # rowc = [big[p]-small[p] for p in range(len(big))].index(1) + 1
            # colc = big[rowc-1]
            return one_step.schutzenberger() + box_evac
        else:
            raise ValueError("must be a semistandard tableau")
            

# def chain_to_tableau(shapes):
#     # converts a list of skewshapes that adds one box at a time into a tableaux
#     # with boxes numbered in order they were created.
#     if shapes[0].boxes != []:
#         fsh = shapes[0].inner_partition
#         shapes = [skewshape(fsh,fsh)] + shapes
#     sh1 = shapes[0]
#     i = 1
#     N = len(sh1.outer_partition)
#     T1 = [[]]*N
#     tabl = tableau(sh1,T1)
#     for sh in shapes[1:]:
#         i += 1
#         N = len(sh.outer_partition)
#         box = sh.boxes[0]
#         T = [[]]*(box[0]-1) + [[i]] + [[]]*(N-box[0])
#         tabl = tabl + tableau(sh,T)
#     return tabl

def chain_to_tableau(shapes):
    # converts a list of paritions to a skew-tableau. The boxes added in the ith step
    # are labelled i.
    if shapes == []:
        raise ValueError("shapes must not be empty")
    inner_partition = shapes[0]
    outer_partition = shapes[-1]
    T = [[]]*len(outer_partition)
    for i in range(1,len(shapes)):
        boxes = skewshape(shapes[i],shapes[i-1]).boxes
        for j in range(1,len(outer_partition)+1):
            T[j-1] = T[j-1] + [i]*len([box for box in boxes if box[0] == j])
    return tableau(skewshape(outer_partition,inner_partition),T)
            

# def tableau_to_chain(T):
#     # converts a skew tableau (must be semistandard) into a chain of shapes
#     if T.is_semistd():
#         shapes = [skewshape(T.shape.inner_partition,T.shape.inner_partition)]
#         for box in T.boxes_ordered:
#             shapes.append(shapes[-1].addbox(box[0],box[1]))
#         return shapes
#     else:
#         raise ValueError("must be a semistandard tableau")

def tableau_to_chain(T):
    # converts a skew tableau (must be semistandard) into a chain of partitions
    if T.is_semisemistd():
        entries = [a for row in T.T for a in row]
        if entries == []:
            return [T.shape.inner_partition,T.shape.outer_partition]
        m = max(entries)
        return [T.shape.inner_partition] + [T.restrict(1,p).shape.outer_partition for p in range(1,m+1)]
    else:
        raise ValueError("must be a semisemistandard tableau")

def YamTab(shape):
    # returns the yamanouchi tableaux of given shape (i.e. 1's in first row, 2's in
    # second row etc)
    T = [[i]*(shape.row_sizes[i-1]) for i in range(1,len(shape.outer_partition)+1)]
    return tableau(shape,T)
        
class tabloid(object):
    def __init__(self,n,T):
        self.T = strip_empties(T)
        self.shape = [len(r) for r in T]
        self.n = n

    def __repr__(self):
        if self.T == []:
            return "empty"
        rep_str = ""
        colwidths = ['.' for i in range(len(self.T[0]))]
        for j in range(len(self.T[0])):
            max_digits = 1
            for entry in self.col(j+1):
                digits = len(str(entry))
                max_digits = max([max_digits,digits])
            colwidths[j] = max_digits
        for row in self.T:
            row_str = ""
            for i in range(len(row)):
                space_buffer = colwidths[i] - len(str(row[i]))
                row_str += str(row[i])+" "*(space_buffer+1)
            rep_str += row_str +"\n"
        return rep_str[:-2]

    def __eq__(self,other):
        return (self.stdrep().T == other.stdrep().T)

    def add_row(self,R):
        #print(R)
        return tabloid(self.n,self.T + [R])

    def rem_row(self):
        return tabloid(self.n,self.T[:-1])

    def col(self,j):
        # returns the jth column
        return [row[j-1] for row in self.T if len(row) >= j]

    def stdrep(self):
        # returns the same tabloid but with rows ordered and reduced mod n (lies in {1,2,...,n}).
        newT = [sorted([aorn(a%self.n,self.n) for a in row]) for row in self.T]
        return tabloid(self.n,newT)

    def skewtab(self):
        # returns the skewtableau defined as follows: take the standard rep of T and
        # slide rows to the right as much as needed so that the skew tableau becomes
        # semistandard/standard.
        std_rows = self.stdrep().T
        row_adj1 = [0]*len(self.T)
        for i in range(1,len(self.T)):
            for p in range(len(self.T[i])+1):
                tests = [(std_rows[i-1][j] < std_rows[i][j+p]) for j in range(len(self.T[i])-p)]
                if all(tests):
                    row_adj1[i-1] = p
                    break
        #row_adj1 = [len([a for a in self.T[i+1] if a < min(self.T[i])]) for i in range(len(self.T)-1)] + [0]
        row_adj = [sum(row_adj1[i:]) for i in range(len(self.T))]
        outer = [row_adj[i]+len(self.T[i]) for i in range(len(self.T))]
        sh = skewshape(outer,row_adj)
        return tableau(sh,std_rows)

    def rect(self,return_evac=False):
        return self.skewtab().rect(return_evac)

    def schutzenberger(self):
        [rec,dec] = self.rect(True)
        return tabloid(self.n,rec.schutzenberger().fwd_slide(dec).T)

# print("starting...")

# T=chain_to_tableau([[],[2],[3,1],[4,3,1],[4,4,1]])
# print(T.schutzenberger())
# T = tabloid(3, [[3],[1],[2]])
# print(T)
# print(T.skewtab())


# 
# sh = skewshape([6,6,5,5,3,2,1],[6,4,4,2,1,1,1])

# T = tableau(sh,[[],[2,3],[3],[1,3,4],[2,3],[4],[]])
# sh2 = skewshape([6,6,6,6,6,2,1],[6,6,5,5,3,2,1])
# S = tableau(sh2,[[],[],[1],[2],[3,4,5],[],[]])

# print(T)
# print(T.fwd_slide_adj(S,True))
#print(S)
#print(S.boxes_ordered)
#print(T.rev_slide_adj(S,True))
# print(T.rev_slide(4,2))

# print(T.rect(True)[0])
# print(T.rect(True)[1])



# print(T.schutzenberger())


# T = tabloid(11,[[6,8,9,10,11],[2,4,5,7],[1,3]])
# print(T)
# print()
# print(T.schutzenberger())

# L = skewshape([4,2,2,2,1],[2,1,1])
# print(L.diagram())
# print(L.outer_addables)
# print(L.outer_removables)
# print(L.inner_addables)
# print(L.inner_removables)

