import matrixOperations
from math import exp, floor

class blackscholes():
    def __init__(self, S, t, K, r, sigma, q):
        self.S = float(S)
        self.t = float(t)
        self.K = float(K)
        self.r = float(r)/100
        self.sigma = float(sigma)
        self.q = float(q)

    #Explicit FDM
    def explicitfdm(self, M, N):
        #1. Compute delta_t = t/N; delta_s = s_max/M, where s_max = 2*K
        delta_t = self.t/N
        s_max = 2 * self.K 
        delta_s = s_max/M

        #2. Compute (Call option) fN,j=max(j*delta_s-K,0), (Put option) fN,j=max(K-j*delta_s), for j=0,1,...,M
        call = [[0 for j in range(M+1)] for x in range(N+1)]
        put = [[0 for j in range(M+1)] for x in range(N+1)]
        for j in range(M + 1):
            call[N][j] = max(j * delta_s - self.K, 0)
            put[N][j] = max(self.K - j * delta_s, 0)

        #3. For i=N-1, N-2,...,1,0, repeat 3.1 and 3.2
        for i in range(N - 1, -1, -1):
            #3.1. Compute vector Fi=A*Fi+1
            for j in range(1, M):
                a = 0.5 * delta_t * (self.sigma * self.sigma * j * j - (self.r - self.q) * j)
                b = 1 - delta_t * (self.sigma * self.sigma * j * j + self.r)
                c = 0.5 * delta_t * (self.sigma * self.sigma * j * j + (self.r - self.q) * j)
                call[i][j] = (a *call[i+1][j-1] + b * call[i+1][j] + c * call[i+1][j+1])
                put[i][j] = (a * put[i+1][j-1] + b * put[i+1][j] + c * put[i+1][j+1])
            #3.2. Compute vector Fi
            call[i][0] = 0
            call[i][M] = s_max - self.K * exp(-self.r * (N-i) * delta_t)
            put[i][0] = self.K * exp(-self.r * (N-i) * delta_t)
            put[i][M] = 0

        #4. Find k, such that k*(delta_s) <= S <= (k+1)*(delta_s), i.e. k=[S/delta_s]
        k = int(floor(self.S / delta_s))

        #5. Option price: V = f0,k + ((f0,k+1 - f0,k) / delta_s) * (S - k*delta_s)
        Vcall = round((call[0][k] + (call[0][k+1] - call[0][k]) / delta_s * (self.S - k * delta_s)),4)
        Vput = round((put[0][k] + put[0][k+1] - put[0][k]) / delta_s * (self.S - k * delta_s),4)
        return [Vcall, Vput]

    #Implicit FDM
    def implicitfdm(self, M, N):  
        delta_t = self.t/N
        s_max = 2 * self.K
        delta_s = s_max / M
        
        callh = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        puth = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        call = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        put = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        for j in range(M + 1):
            call[j][N] = max(j * delta_s - self.K, 0)
            put[j][N] = max(self.K - j * delta_s, 0)
        
        A = [[0 for _ in range(M + 1)] for _ in range(M + 1)]
        A[0][0] = 1
        A[M][M] = 1
        for j in range(1, M):
            A[j][j - 1] = 0.5 * delta_t * ((self.r - self.q) * j - self.sigma * self.sigma * j * j)
            A[j][j] = 1 + delta_t * (self.sigma * self.sigma * j * j + self.r - self.q)
            A[j][j + 1] = -0.5 * delta_t * (self.sigma * self.sigma * j * j + (self.r - self.q) * j)
        A_inv = matrixOperations.inverse(A)
        for i in range(N - 1, -1, -1):
            # 3.1 0 for j=0, same for j=1toM-1, exponential discounting for highest payout
            for j in range(M+1):
                callh[j][i + 1] = call[j][i + 1]
                puth[j][i + 1] = put[j][i + 1]
            callh[0][i + 1] = 0
            callh[M][i + 1] = s_max - self.K * exp(-self.r * (N - i) * delta_t)
            puth[M][i + 1] = 0
            puth[0][i + 1] = self.K * exp(-self.r * (N - i) * delta_t)
            # 3.2 store vectors from temp cols
            callh_cv = [[callh[j][i + 1]] for j in range(M + 1)]
            puth_cv = [[puth[j][i + 1]] for j in range(M + 1)]
            call_cv = matrixOperations.multiply(A_inv, callh_cv)
            put_cv = matrixOperations.multiply(A_inv, puth_cv)
            for j in range(M + 1):
                call[j][i] = call_cv[j][0]
                put[j][i] = put_cv[j][0]
        
        k = int(floor(self.S / delta_s))

        Vcall = round((call[k][0] + (call[k + 1][0] - call[k][0]) / delta_s * (self.S - k * delta_s)),4)
        Vput = round((put[k][0] + (put[k + 1][0] - put[k][0]) / delta_s * (self.S - k * delta_s)),4)
        return [Vcall, Vput]

