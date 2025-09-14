def gcd(x,y)
    return y==0?x:gcd(y,x%y)
end

def lcm(x,y)
    return x*y/gcd(x,y)
end

def modpow(x,n,mod)
    res = 1
    while n>0 do
        res = res*x%mod if n&1==1
        x = x*x%mod
        n>>=1
    end
    return res
end

def sqrt_floor(n)
    s=1; 
    while (s*s>n||n>=(s+1)*(s+1)) do 
        s=(n/s+s)/2
    end 
    return s
end

def g_dist(x1,y1,x2,y2)
    return Math.sqrt((x1-x2)**2 + (y1-y2)**2)
end

def m_dist(x1,y1,x2,y2)
    return (x1-x2).abs + (y1-y2).abs
end

class PriorityQueue
    @heap
    @size
    @reverse
    def initialize(str)
        @heap = Array.new()
        @size = 0
    end

    def push(x)
        x = x
        i = @size
        @size+=1
        @heap.append(x)
        while i>0 do
            par = (i-1)/2
            break if(@heap[par]<=x)
            @heap[i] = @heap[par]
            i = par
        end

        @heap[i] = x
    end

    def pop()
        ret = @heap[0]
        @size -= 1
        x = @heap[@size]
        i = 0
        while 2*i + 1<@size do
            a = 2*i+1
            b = 2*i+2
            a = b if b<@size and @heap[b] < @heap[a]
            break if @heap[a] >= x

            @heap[i] = @heap[a]
            i = a
        end

        @heap[i] = x
        @heap.pop

        return ret
    end

    def top()
        return @heap[0]
    end

    def size()
        return @size
    end

    def empty()
        return @size == 0
    end
end

class SegmentTree
    @dat
    @sz
    @mono
    def initialize(n,mono)
        @sz = 1
        while @sz < n do
            @sz*=2
        end
        @mono = mono
        @dat = Array.new(2*@sz-1){@mono}
    end

    def f(a,b)
        return gcd(a,b)
    end

    def g(a,b)
        return a+b
    end

    def setval(i,x)
        i += @sz - 1
        @dat[i] = x
        while i>0 do
            i = (i-1)/2
            @dat[i] = f(@dat[2*i+1],@dat[2*i+2])
        end
    end
    
    def update(i,x)
        i += @sz - 1
        @dat[i] = g(@dat[i],x)
        while i>0 do
            i = (i-1)/2
            @dat[i] = f(@dat[2*i+1],@dat[2*i+2])
        end
    end

    def query(l,r)
        l += @sz-1
        r += @sz-1
        vl = @mono.clone
        vr = @mono.clone
        while l<=r do
            vl = f(vl,@dat[l]) if l&1==0
            vr = f(@dat[r],vr) if r&1==1
            l>>=1
            r = (r>>1)-1
        end
        return f(vl,vr)
    end
end

class UnionFind
    @par
    def initialize(n)
        @par = Array.new(n,-1)
    end

    def root(x)
        return x if @par[x]<0
        return @par[x] = root(@par[x])
    end

    def size(x)
        return -@par[root(x)]
    end

    def issame(a,b)
        return root(a)==root(b)
    end

    def connect(a,b)
        a = root(a)
        b = root(b)
        return false if a==b
        if size(a)<size(b) then a,b = b,a end
        @par[a] += @par[b]
        @par[b] = a
        return true
    end
end