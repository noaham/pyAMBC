from itertools import combinations
from math import ceil
from tableaux import *

# This is a dummy variable
x = 'x'

# takes a number and reduces it to it's rep mod n, in 1,2,...,n
def modn(i,n):
    if i%n == 0:
        return n
    else:
        return i%n

# this class defines a ball. Relations are given by NW ordering. So (0,0) <= (1,1).
class ball(object):
    def __init__(self,ii,jj):
        self.i = int(ii)
        self.j = int(jj)

    def __eq__(self, other):
        return ((self.i, self.j) == (other.i, other.j))

    def __ne__(self, other):
        return ((self.i, self.j) != (other.i, other.j))

    def __le__(self, other):
        return ((self.i <= other.i) and (self.j <= other.j))

    def __ge__(self, other):
        return ((self.i >= other.i) and (self.j >= other.j))

    def __lt__(self, other):
        return  ((self <= other) and (self != other))

    def __gt__(self, other):
        return ((self >= other) and (self != other))

    def __repr__(self):
        return "(%s %s)" % (self.i, self.j)

    def ht(self):
        return self.j-self.i

    def nrow(self,n):
        # divide the matrix into nxn grid. Row 1 contains the rows 1..n.
        # returns the row the ball is in (i.e. returns 1 if ball.i = 1..n)
        return int((self.i - modn(self.i,n))/n + 1)

    def translate(self,t):
        # translates by (t,t)
        return ball(self.i+t,self.j+t)

    def trans_row1(self,n):
        # finds translate in the first row. We need to make the nrow 1, so if it is
        # in nrow k, then we need to shift by -(k-1)n
        return self.translate(-(self.nrow(n)-1)*n)

# This sorts in the NE order (useful for zig zags)
def sort_NE(balls):
    flipped = [ball(b.i,-b.j) for b in balls]
    return [ball(b.i,-b.j) for b in sorted(flipped)]
    

# We use the following meathod to encode partial permutations:
#   [4,x,x,1,x,8]
# in the usual window notation. In the above case p(1)=4, p(2) is undefined and p(7)=4+6=10

# This converts the window notation to a set of balls
def window_to_balls(window):
    balls = []
    for i in range(len(window)):
        if window[i] != 'x':
            balls.append(ball(i+1,window[i]))
    return balls

def balls_to_window(balls,n):
    window = [x for i in range(n)]
    for b in balls:
        if b.i in range(1,n+1):
            window[b.i-1] = b.j
        else:
            t = ((b.i % n) - b.i)/n
            if b.i%n == 0:
                bnew = ball(b.i+(t+1)*n,b.j+(t+1)*n)
            else:
                bnew = ball(b.i+t*n,b.j+t*n)
            window[bnew.i-1] = bnew.j
    return window

# This determines whether a set of balls is a chain in the NW ordering
def is_chain(balls):
    if len(balls)==0:
        return True
    a = balls[0]
    for b in balls[1:]:
        if (not ((a <= b) or (b <= a))) or a == b :
            return False
    return is_chain(balls[1:])

class stream(object):
    # defines a class representing a stream. Given a set of balls, all n-translates
    # should form a chain.
    def __init__(self,balls,n):
        self.n = n
        first_row = []
        self.balls = sorted([ball.trans_row1(self.n) for ball in balls])
        self.k = len(self.balls)
        if is_chain(self.balls):
            if self.k != 0:
                M = max(self.balls)
                m = min(self.balls)
                if not ((M.translate(-n) <= m) and (M <= m.translate(n))):
                    raise ValueError("balls must form a chain")
        else:
            raise ValueError("balls must form a chain")
        [self.A,self.B,self.r] = self.data()

    def __eq__(self,other):
        if self.n == other.n:
            return [self.A,self.B,self.r] == [other.A,other.B,other.r]
        else:
            raise ValueError("streams must have same n to be comparable")

    def __repr__(self):
        return str(self.balls)

    def __contains__(self,ball):
        # allows us to write "ball in stream" and get a True/False result
        return ball.trans_row1(self.n) in self.balls

    def __getitem__(self,i):
        # allows us to write stream[i] and get the ith ball in the stream.
        if self.k == 0:
            raise ValueError("string must be nonempty to get element")
        r = modn(i,self.k)
        row1_ball = self.balls[r-1]
        t = (i-r)/self.k
        return row1_ball.translate((t)*self.n)
    
    def ball(self,i):
        # returns the ith ball in the stream (first ball is the one
        # just south east of (0,0) and number increases southeast)
        # we have that i = r + tk for some r and t. t will be one less than the
        # row we need to look in and it will be the rth ball in that row.
        if self.k == 0:
            raise ValueError("this function needs k to be nonzero")
        r = modn(i,self.k)
        row1_ball = self.balls[r-1]
        t = (i-r)/self.k
        return row1_ball.translate((t)*self.n)

    def rowi_balls(self,i):
        # returns the balls in row i
        return [ball.translate((i-1)*self.n) for ball in self.balls]

    def data(self):
        # returns A,B,r where A are the row residues, B the column residues
        # and r is the altitude
        if self.k == 0:
            return [[],[],None]
        A = sorted([modn(ball.i,self.n) for ball in self.balls])
        B = sorted([modn(ball.j,self.n) for ball in self.balls])
        i = self.ball(1).i
        j = self.ball(1).j
        # how many columns over is j?
        mv = int((j - modn(j,self.n))/self.n)
        r = B.index(modn(j,self.n)) + mv*self.k
        return [A,B,r]

    def nearest_NW(self,ball):
        # returns the ball b in stream that is NW of ball and maximumlly so
        if self[1] <= ball:
            i=0
            while self[1+i*self.k] <= ball:
                i += 1
            checks = self.rowi_balls(i-1) + self.rowi_balls(i)
            return max([b for b in checks if b <= ball])
        else:
            i=0
            while not self[1-i*self.k] <= ball:
                i += 1
            checks = self.rowi_balls(-i+1)
            return max([b for b in checks if b <= ball])
        
    
    def d(self,ball):
        # returns the numbering of the stream, i.e. d(balls[0]) = 1 increasing SE
        # however if the ball is not in the stream then it returns d(b) where b is the
        # maximum ball in the stream that occurs NW of ball.
        if ball in self:
            row = ball.nrow(self.n)
            return self.balls.index(ball.trans_row1(self.n))+1+(row-1)*self.k
        else:
            return self.d(self.nearest_NW(ball))

def st_r(A,B,r,n):
    # returns the stream str_r(A,B)
    Amod = sorted([modn(a,n) for a in A])
    Bmod = sorted([modn(b,n) for b in B])
    k = len(Amod)
    if k == len(set(Amod)) and k == len(set(Bmod)):
        balls = []
        for i in range(1,k+1):
            j = modn(i+r,k)
            # i+r = j+tk
            t = int((i+r-j)/k)
            balls.append(ball(Amod[i-1],Bmod[j-1]+t*n))
        return stream(balls,n)
    else:
        raise ValueError("A and B must be equal size")

def is_stream(balls,n):
    if is_chain(balls):
        M = max(balls)
        m = min(balls)
        return ((ball(M.i - n,M.j - n) <= m) and (M <= ball(m.i+n,m.j+n)))
    return False
    
# this function calculates all the balls that are reachable in the NW order from a given ball b. In fact,
# it calculates the first translate of each ball that is reachable.
def reachable_balls(balls,n,b):
    reaches = []
    for a in balls:
        if a.trans_row1(n) != b.trans_row1(n):
            t = 0
            if a < b:
                while a.translate(t*n) < b:
                    t += 1
                reaches.append(a.translate((t-1)*n))
            else:
                while not (a.translate(-t*n) < b):
                    t += 1
                reaches.append(a.translate(-t*n))
    return reaches

# checks if b is an n translate of a stream
def is_in_stream(n,stream,b):
    stream_row = int(ceil((stream[0].i)/float(n)))
    b_row = int(ceil((b.i)/float(n)))
    t = b_row - stream_row
    return ball(b.i-t*n,b.j-t*n) in stream
        

# Here we calculate all paths from a given ball b to a given stream (in practice always a channel)
# restrictions on the paths are that they have at most #balls elements, don't repeat translation classes
# and do not skip intermediate balls. This reduces it to a finite problem. In addition we can rule out paths
# that intersect the stream.
def paths_from_b_to_stream_length(balls,n,b,stream,length):
    paths = []
    reaches = reachable_balls(balls,n,b)
    if length == 0:
        if is_in_stream(n,stream,b):
            return [[b]]
        else:
            return []
    elif length > 0:
        for a in reaches:
            inbetweens = [((a < c) and (c < b)) for c in reaches]
            is_neighbour = not any(inbetweens)
            not_in_stream = not is_in_stream(n,stream,a)
            if length == 1:
                if (is_neighbour and (not not_in_stream)):
                    paths.append([b,a])
                    
            else:
                if is_neighbour and not_in_stream:
                    paths_at_a = paths_from_b_to_stream_length(balls,n,a,stream,length-1)
                    paths = paths + [[b]+path for path in paths_at_a]
    return paths

            

class partialperm(object):
    def __init__(self, windoww):
        self.window = windoww
        self.n = len(self.window)
        self.balls = window_to_balls(self.window)
        self.num_balls = len(self.balls)
        self.streams_data = self.streams_width_channels()
        self.streams = self.streams_data[0]
        self.width = self.streams_data[1]
        self.channels = self.streams_data[2]
        if self.balls == []:
            min_ht = 0
            min_chs = []
        else:
            min_ht = min([sum([b.ht() for b in ch]) for ch in self.channels])
            min_chs = [ch for ch in self.channels if sum([b.ht() for b in ch])==min_ht]
        if len(min_chs)==1:
            self.SW_channel = sorted(min_chs[0])
        else:
            self.SW_channel = []
        self.numbering = [self.ch_numbering(b) for b in self.balls]
        self.zigzags = self.zig_zag_constr()
        self.back_posts = [x] + [ball(zz[0].i,zz[-1].j) for zz in self.zigzags[1:]]
        self.inner_corners_w = balls_to_window(self.inner_corners(),self.n)
        self.st_r = stream(self.back_posts[1:],self.n)
        
        self.Prow = sorted([modn(b.j,self.n) for b in self.back_posts[1:]])
        self.Qrow = sorted([modn(b.i,self.n) for b in self.back_posts[1:]])
        self.rho_row = int(sum([ceil(b.j/float(self.n))-ceil(b.i/float(self.n)) for b in self.back_posts[1:]]))

    def __repr__(self):
        repstr = '[ '
        for i in self.window:
            repstr += str(i)+' '
        return repstr + ']'
        
    def streams_width_channels(self):
        is_width_def = False
        streams = []
        if len(self.balls) == 0:
            width = 0
            channels = []
        for d in range(1,len(self.balls)+1)[::-1]:
            subs = [list(L) for L in combinations(self.balls,d)]
            for subset in subs:
                if is_stream(subset,self.n):
                    streams.append(sorted(subset))
            if ((not is_width_def) and (len(streams) > 0)):
                channels = [S for S in streams]
                width = len(streams[0])
                is_width_def = True
        return [streams,width,channels]

    def ch_numbering(self,b):
        if b in self.balls:
            if b in self.SW_channel:
                return self.SW_channel.index(b)+1
            else:
                paths = []
                for i in range(1,len(self.balls)):
                    paths = paths + paths_from_b_to_stream_length(self.balls,self.n,b,self.SW_channel,i)
                max_worth = max([self.ch_numbering(path[-1])+len(path)-1 for path in paths])
                return max_worth
        else:
            b_row = int(ceil((b.i)/float(self.n)))
            b_trans = ball(b.i-(b_row-1)*self.n,b.j-(b_row-1)*self.n)
            return self.ch_numbering(b_trans)+(self.width)*(b_row-1)

    def zig_zag_constr(self):
        zigzags = [x for i in range(self.width+1)]
        for z in range(1,self.width+1):
            zag = []
            for i in range(len(self.balls)):
                if (self.numbering[i] % self.width) == (z % self.width):
                    b = self.balls[i]
                    d = self.numbering[i]
                    t = (z-d)/self.width
                    trans_b = ball(b.i+t*self.n,b.j+t*self.n)
                    if self.ch_numbering(trans_b) == z:
                        zag.append(trans_b)
            zigzags[z] = sort_NE(zag)
        return zigzags

    def rowi_balls(self,i):
        # returns the balls in row i
        return [ball.translate((i-1)*self.n) for ball in self.balls]
    
    def inner_corners(self):
        inners = []
        for zz in self.zigzags[1:]:
            for i in range(1,len(zz)):
                topb = zz[i-1]
                botb = zz[i]
                inners.append(ball(botb.i,topb.j))
        return inners

    def is_monotone(self,d,example=False):
        # input is a numbering d of the balls and returns whether or not it
        # is monotone
        for b in self.balls:
            for a in [ball for ball in self.balls if ball != b]:
                ar = reachable_balls([a],self.n,b)[0]
                if d(ar) == d(b):
                    return b if example else False
        return True
        
    def bk_stream_numbering(self,stream):
        # given a stream (rows determined by A and cols by B) of altitude r, determine the
        # backward numbering of the parital perm
        nums = [stream.d(ball) for ball in self.balls]
        k = stream.k
        d = lambda ball : nums[self.balls.index(ball.trans_row1(self.n))] + k*(ball.nrow(self.n)-1)
        while not self.is_monotone(d):
            exb = self.is_monotone(d,True)
            while [a for a in reachable_balls(self.balls,self.n,exb) if d(a)==d(exb)] != []:
                exb = [a for a in reachable_balls(self.balls,self.n,exb) if d(a)==d(exb)][0]
            i = self.balls.index(exb.trans_row1(self.n))
            nums[i] -= 1
        return nums

    def bk_w(self,stream):
        # using the backward stream numbering, create zig zags with outer corners on the balls
        # in self and back posts being the stream. The inner corners define a partial perm
        # which is returned.
        zigzags = []
        nums = self.bk_stream_numbering(stream)
        k = stream.k
        d = lambda ball : nums[self.balls.index(ball.trans_row1(self.n))] + k*(ball.nrow(self.n)-1)
        for b in stream.balls:
            outers = [a.translate(int((stream.d(b)-d(a))*self.n/k)) for a in self.balls if stream.d(b)%k == d(a)%k]
            outers = sort_NE(outers)
            zz = []
            i = b.i
            for a in outers:
                j = a.j
                zz.append(ball(i,j))
                i = a.i
            zz.append(ball(i,b.j))
            zigzags.append(zz)
        bk_w_balls = [b.trans_row1(self.n) for zz in zigzags for b in zz]
        return partialperm(balls_to_window(bk_w_balls,self.n))
            
        
        # if len(A) == len(B):
        #     k = len(A)
        #     Amod = ordered([modn(a,self.n) for a in A])
        #     Bmod = ordered([modn(b,self.n) for b in B])
        #     if k != len(set(Amod)) or k != len(set(Bmod)):
        #         raise ValueError("A and B must contain at most one representative of each mod n equivalence class")
        #     stream_balls0 = [ball(Amod[i],Bmod[i]) for i in range(k)]
        #     stream_ballsr = 0
        # else:
        #     raise ValueError("A and B must be the same size")
        


def AMBC(window):
    perm = partialperm(window)
    P = tabloid(perm.n,[])
    Q = tabloid(perm.n,[])
    rho = [] # tabloid(perm.n,[])
    while perm.num_balls > 0:
        [Prow,Qrow,rho_row] = perm.st_r.data()
        P = P.add_row(perm.Prow)
        Q = Q.add_row(perm.Qrow)
        rho.append(perm.rho_row) # rho.add_row([perm.rho_row])
        perm = partialperm(perm.inner_corners_w)
    return [P,Q,rho]

def invAMBC(P,Q,rho):
    P = P.stdrep()
    Q = Q.stdrep()
    n = P.n
    window = ['x']*n
    perm = partialperm(window)
    while len(P.T) > 0:
        B = P.T[-1]
        A = Q.T[-1]
        r = rho[-1]
        bk_str = st_r(A,B,r,n)
        perm = perm.bk_w(bk_str)
        P = P.rem_row()
        Q = Q.rem_row()
        rho = rho[:-1]
    return perm.window

def RSK(window):
    # window must be a word
    if window == []:
        sh = skewshape([])
        return (tableau(sh,[]),tableau(sh,[]))
    (Pprev,Qprev) = RSK(window[:-1])
    P = Pprev*tableau(skewshape([1]),[[window[-1]]])
    box = skewshape(P.shape.outer_partition,Pprev.shape.outer_partition)
    (i,j) = box.boxes[0]
    Q = Qprev.addbox(i,j,len(window))
    return (P,Q)
    
    # w = window
    # n = len(w)
    # outer = [n-i for i in range(n)]
    # inner = [a-1 for a in outer]
    # sh = skewshape(outer,inner)
    # word_tabl = tableau(sh,[[a] for a in w][::-1])
    # P = word_tabl.rect()
    # winv = [w.index(i)+1 for i in range(1,n+1)]
    # word_tabl = tableau(sh,[[a] for a in winv][::-1])
    # Q = word_tabl.rect()
    # return [P,Q]

def invRSK(P,Q):
    # inverse of the rsk algorithm. We use the reverse jdt procedure
    if P.shape == Q.shape and P.is_semistd() and Q.is_std():
        sh = P.shape
        if sh == skewshape([]):
            return []
        # locate the largest entry in Q
        (i,j) = Q.boxes_ordered[-1]
        Qprev = Q.rembox(i,j,'outer')
        slides = [(len(P.col(l))+1,l) for l in range(1,sh.outer_partition[0]+1) if l != j]
        slides.append((1,sh.outer_partition[0]+1))
        tempsh = skewshape(sh.outer_partition,sh.outer_partition)
        T = tableau(tempsh,[[]]*len(sh.outer_partition))
        c = 1
        for (k,l) in slides:
            T = T.addbox(k,l,c)
            c += 1
        Sep = P.fwd_slide(T)
        letter = Sep.T[0][-1]
        dj = Sep.shape.outer_partition[0]
        Pprev = Sep.rembox(1,dj).rect()
        return (invRSK(Pprev,Qprev) + [letter])
    else:
        raise ValueError("P and Q must be the same shape and P must be semistandard while Q is standard")

#print("starting...")
# w = [49, 22, 5, 23, 52, 10, 19, 24]
# perm = partialperm(w)
# print(AMBC(w)[0])
# print()
# print(AMBC(w)[1])
# print()
# print(AMBC(w)[2])

# sh = skewshape([2,1])
# A = tableau(sh,[[1,2],[3]])
# B = tableau(sh,[[1,3],[2]])
# print(RSK(invRSK(A,B)))

# s1 = stream([ball(1,-1),ball(3,2)],4)
# s2 = stream([ball(1,3),ball(3,6)],4)

#print(s2.nearest_NW(ball(1,3)))
#s = st_r([2,3,5,7],[1,4,5,7],2,7)
#print(s0)
# p = partialperm([3,'x','x',9,'x',6,'x'])
# for b in p.balls:
#     print(b,p.bk_stream_numbering(s)[p.balls.index(b)])
# print(p.bk_w(s))
# print(invAMBC(tabloid(7,[[2,3,5,7],[1,4],[6]]),tabloid(7,[[1,4,5,7],[3,6],[2]]),[2,0,1]))
# print(p.balls)
# print(p.bk_stream_numbering(s0))
#print(p.balls.index(b.trans_row1(p.n)))
# for b in p.rowi_balls(2):
#     print(b)
#     print(s0.d(b))
#     print(p.bk_stream_numbering(s0)(b))

# print(s0.nearest_NW(ball(2,6)))
# print(p.width)
