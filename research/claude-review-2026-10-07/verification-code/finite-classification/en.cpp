// Independent verification of terminal-element classification in E6, E7, E8.
// Written from scratch for the referee report. Exact integer arithmetic.
// Numbering: chain 0..n-2, node n-1 attached to node 2.
// Element = matrix of w(alpha_j) in simple-root coordinates, column j = w(alpha_j).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <functional>
using namespace std;

static int N;                 // rank
static int C[8][8];           // Cartan matrix
static bool ADJ[8][8];

struct Mat {
    int8_t a[8][8]; // a[col][row]: coefficient of alpha_row in w(alpha_col)
    bool operator<(const Mat& o) const { return memcmp(a, o.a, sizeof(a)) < 0; }
    bool operator==(const Mat& o) const { return memcmp(a, o.a, sizeof(a)) == 0; }
};
struct MatHash { size_t operator()(const Mat& m) const {
    uint64_t h = 1469598103934665603ULL; const uint8_t* p=(const uint8_t*)m.a;
    for (size_t i=0;i<sizeof(m.a);i++){ h ^= p[i]; h *= 1099511628211ULL; } return h; } };

static void setup(int n) {
    N = n;
    memset(C,0,sizeof(C)); memset(ADJ,0,sizeof(ADJ));
    for (int i=0;i<n;i++) C[i][i]=2;
    auto edge=[&](int i,int j){ ADJ[i][j]=ADJ[j][i]=true; C[i][j]=C[j][i]=-1; };
    for (int i=0;i+1<=n-2;i++) edge(i,i+1);
    edge(2,n-1);
}
static Mat identity() { Mat m; memset(m.a,0,sizeof(m.a)); for(int i=0;i<N;i++) m.a[i][i]=1; return m; }
// w -> w s_i : column i negated, columns j~i get + old column i
static void rightMul(Mat& m, int i) {
    int8_t old[8]; memcpy(old, m.a[i], 8);
    for (int k=0;k<N;k++) m.a[i][k] = -old[k];
    for (int j=0;j<N;j++) if (ADJ[i][j]) for (int k=0;k<N;k++) m.a[j][k] += old[k];
}
// w -> s_i w : each column v -> v - B(v,alpha_i) alpha_i
static void leftMul(Mat& m, int i) {
    for (int j=0;j<N;j++) {
        int p=0; for (int k=0;k<N;k++) p += m.a[j][k]*C[k][i];
        m.a[j][i] -= p;
    }
}
static Mat mul(const Mat& A, const Mat& B) { // (A*B)(alpha_j) = A(B(alpha_j))
    Mat R; memset(R.a,0,sizeof(R.a));
    for (int j=0;j<N;j++) for (int k=0;k<N;k++) { int bkj = B.a[j][k]; if(!bkj) continue;
        for (int r=0;r<N;r++) R.a[j][r] += bkj * A.a[k][r]; }
    return R;
}
static bool positive(const int8_t* v) { bool nz=false; for(int k=0;k<N;k++){ if(v[k]<0) return false; if(v[k]>0) nz=true;} return nz; }
static bool negative(const int8_t* v) { bool nz=false; for(int k=0;k<N;k++){ if(v[k]>0) return false; if(v[k]<0) nz=true;} return nz; }
static bool isRightDescent(const Mat& m, int i) { return negative(m.a[i]); }
// right-terminal test (paper eq. terminaltest): s in R(w) => w(alpha_s+alpha_t) > 0 for all t~s
static bool rightTerminal(const Mat& m) {
    for (int s=0;s<N;s++) if (isRightDescent(m,s)) {
        for (int t=0;t<N;t++) if (ADJ[s][t]) {
            int8_t v[8]; for(int k=0;k<N;k++) v[k]=m.a[s][k]+m.a[t][k];
            if (!positive(v)) return false;
        }
    }
    return true;
}
// alternative right-terminal test: no s in R(w), t~s with t in R(ws)
static bool rightTerminal2(const Mat& m) {
    for (int s=0;s<N;s++) if (isRightDescent(m,s)) {
        Mat ws=m; rightMul(ws,s);
        for (int t=0;t<N;t++) if (ADJ[s][t] && isRightDescent(ws,t)) return false;
    }
    return true;
}
// exact inverse via w^{-1} = C^{-1} w^T C  (B-orthogonality), computed with adj(C)/det(C)
static long long detC; static long long adjC[8][8];
static long long det(vector<vector<long long>> M) {
    int n=M.size(); long long d=1;
    // fraction-free Bareiss
    long long prev=1;
    for (int k=0;k<n-1;k++) {
        if (M[k][k]==0) { int sw=-1; for(int r=k+1;r<n;r++) if(M[r][k]!=0){sw=r;break;} if(sw<0) return 0; swap(M[k],M[sw]); d=-d; }
        for (int i=k+1;i<n;i++) for (int j=k+1;j<n;j++) M[i][j] = (M[i][j]*M[k][k]-M[i][k]*M[k][j])/prev;
        prev=M[k][k];
    }
    return d*M[n-1][n-1];
}
static void setupInverse() {
    vector<vector<long long>> M(N, vector<long long>(N));
    for(int i=0;i<N;i++) for(int j=0;j<N;j++) M[i][j]=C[i][j];
    detC = det(M);
    for (int i=0;i<N;i++) for (int j=0;j<N;j++) {
        vector<vector<long long>> S; for(int r=0;r<N;r++){ if(r==j) continue; vector<long long> row; for(int c=0;c<N;c++){ if(c==i) continue; row.push_back(C[r][c]);} S.push_back(row);}
        long long cof = (N==1)?1:det(S); if ((i+j)&1) cof=-cof;
        adjC[i][j] = cof; // adj = transpose of cofactor matrix; C symmetric so fine
    }
}
static Mat inverse(const Mat& m) {
    // w as matrix W[row][col] = m.a[col][row]. w^{-1} = C^{-1} W^T C
    long long T[8][8]; // W^T C : T[r][c] = sum_k W[k][r] C[k][c] = sum_k m.a[r][k] C[k][c]
    for (int r=0;r<N;r++) for (int c=0;c<N;c++){ long long s=0; for(int k=0;k<N;k++) s += (long long)m.a[r][k]*C[k][c]; T[r][c]=s; }
    Mat R;
    for (int r=0;r<N;r++) for (int c=0;c<N;c++){ long long s=0; for(int k=0;k<N;k++) s += adjC[r][k]*T[k][c];
        if (s % detC != 0) { fprintf(stderr,"inverse not integral!\n"); exit(1); }
        long long v = s/detC; if (v<-127||v>127){fprintf(stderr,"overflow\n"); exit(1);} R.a[c][r]=(int8_t)v; }
    return R;
}
static bool isCommutingProduct(const Mat& m) {
    // w is a product of commuting generators iff R(w) is commuting and prod_{s in R(w)} s == w
    Mat t = identity(); vector<int> R;
    for (int s=0;s<N;s++) if (isRightDescent(m,s)) R.push_back(s);
    for (size_t i=0;i<R.size();i++) for (size_t j=i+1;j<R.size();j++) if (ADJ[R[i]][R[j]]) return false;
    for (int s: R) rightMul(t,s);
    return t==m;
}
static string word(Mat m) { // reduced word by stripping right descents (smallest label first)
    string w;
    while (true) { int s=-1; for(int i=0;i<N;i++) if(isRightDescent(m,i)){s=i;break;} if(s<0) break; w.push_back('0'+s); rightMul(m,s); }
    reverse(w.begin(), w.end()); return w;
}
static int length(const Mat& m) { return (int)word(m).size(); }
static string descents(const Mat& m, const Mat& inv) {
    string L,R; for(int i=0;i<N;i++){ if(isRightDescent(inv,i)) L.push_back('0'+i); if(isRightDescent(m,i)) R.push_back('0'+i);} return "L="+L+" R="+R;
}
static string support(const string& w){ set<char> s(w.begin(),w.end()); return string(s.begin(),s.end()); }

// ---------- full enumeration by level BFS (E6, E7) ----------
static void fullEnumerate(bool countFC) {
    vector<Mat> level; level.push_back(identity());
    vector<long long> poincare;
    long long total=0, rightTerm=0, twoSided=0, commT=0, noncommT=0;
    vector<string> noncommWords;
    // FC via recurrence: FC_k set
    unordered_set<Mat,MatHash> fcPrev; if (countFC) fcPrev.insert(identity());
    long long fcTotal = countFC?1:0;
    int len=0;
    while (!level.empty()) {
        poincare.push_back(level.size()); total += level.size();
        for (const Mat& m : level) {
            bool rt = rightTerminal(m);
            if (rt != rightTerminal2(m)) { fprintf(stderr,"terminal tests disagree\n"); exit(1); }
            if (rt) { rightTerm++; Mat inv = inverse(m); if (rightTerminal(inv)) { twoSided++;
                if (isCommutingProduct(m)) commT++; else { noncommT++; noncommWords.push_back(word(m)+"  len="+to_string(len)+" "+descents(m,inv)); } } }
        }
        // next level
        vector<Mat> next; next.reserve(level.size()*N);
        for (const Mat& m : level) for (int i=0;i<N;i++) if (!isRightDescent(m,i)) { Mat t=m; rightMul(t,i); next.push_back(t); }
        sort(next.begin(), next.end()); next.erase(unique(next.begin(), next.end()), next.end());
        if (countFC) {
            unordered_set<Mat,MatHash> fcNext;
            for (const Mat& u : next) {
                vector<int> R; for(int s=0;s<N;s++) if(isRightDescent(u,s)) R.push_back(s);
                bool ok=true;
                for (size_t i=0;i<R.size()&&ok;i++) for (size_t j=i+1;j<R.size();j++) if (ADJ[R[i]][R[j]]) {ok=false;break;}
                for (size_t i=0;i<R.size()&&ok;i++) { Mat t=u; rightMul(t,R[i]); if (!fcPrev.count(t)) ok=false; }
                if (ok) fcNext.insert(u);
            }
            fcTotal += fcNext.size(); fcPrev.swap(fcNext);
        }
        level.swap(next); len++;
    }
    printf("E%d: |W|=%lld maxlen=%d\n", N, total, (int)poincare.size()-1);
    printf("Poincare:"); for (auto c: poincare) printf(" %lld", c); printf("\n");
    printf("right-terminals=%lld two-sided terminals=%lld (commuting %lld, noncommuting %lld)\n", rightTerm, twoSided, commT, noncommT);
    for (auto& w: noncommWords) printf("  noncommuting terminal: %s\n", w.c_str());
    if (countFC) printf("FC count=%lld\n", fcTotal);
}

// ---------- FC enumeration only (closure under recurrence) ----------
static long long fcEnumerate() {
    unordered_set<Mat,MatHash> prev; prev.insert(identity()); long long total=1; int maxlen=0; int len=0;
    while (!prev.empty()) {
        unordered_set<Mat,MatHash> cand;
        for (const Mat& m: prev) for (int i=0;i<N;i++) if (!isRightDescent(m,i)) { Mat t=m; rightMul(t,i); cand.insert(t); }
        unordered_set<Mat,MatHash> nxt;
        for (const Mat& u: cand) {
            vector<int> R; for(int s=0;s<N;s++) if(isRightDescent(u,s)) R.push_back(s);
            bool ok=true;
            for (size_t i=0;i<R.size()&&ok;i++) for (size_t j=i+1;j<R.size();j++) if (ADJ[R[i]][R[j]]) {ok=false;break;}
            for (size_t i=0;i<R.size()&&ok;i++) { Mat t=u; rightMul(t,R[i]); if (!prev.count(t)) ok=false; }
            if (ok) nxt.insert(u);
        }
        total += nxt.size(); if(!nxt.empty()) maxlen=len+1; prev.swap(nxt); len++;
    }
    printf("E%d FC count (closure)=%lld, max FC length=%d\n", N, total, maxlen);
    return total;
}

// ---------- parabolic pruning for E8 ----------
// Build chain by adding nodes in given order. T_J = right-terminals of W_J (ambient matrices).
static void pruning(const vector<int>& order) {
    vector<Mat> T; T.push_back(identity());
    vector<int> J;
    long long totalCandidates=0;
    for (size_t step=0; step<order.size(); step++) {
        int k = order[step]; vector<int> Jp = J; Jp.push_back(k);
        // minimal left coset reps a in W_{J'} with a(alpha_j)>0 for j in J: orbit of omega_k under W_{J'} (weight coords)
        struct Node { array<int,8> lam; Mat a; };
        vector<Node> reps; map<array<int,8>, int> seen;
        Node start; start.lam.fill(0); start.lam[k]=1; start.a=identity();
        reps.push_back(start); seen[start.lam]=0;
        for (size_t idx=0; idx<reps.size(); idx++) {
            Node cur = reps[idx];
            for (int i : Jp) if (cur.lam[i] > 0) {
                Node nx = cur; // s_i lam : lam_j -= lam_i C[i][j]
                for (int j=0;j<N;j++) nx.lam[j] = cur.lam[j] - cur.lam[i]*C[i][j];
                leftMul(nx.a, i);
                if (!seen.count(nx.lam)) { seen[nx.lam]=reps.size(); reps.push_back(nx); }
            }
        }
        // sanity: each rep has no right descent in J, reps distinct as matrices
        set<Mat> distinct;
        for (auto& r: reps) { for (int j: J) if (isRightDescent(r.a,j)) { fprintf(stderr,"rep has J-descent\n"); exit(1);} distinct.insert(r.a); }
        if (distinct.size()!=reps.size()) { fprintf(stderr,"reps not distinct\n"); exit(1); }
        // candidates
        vector<Mat> Tn;
        for (auto& r: reps) for (const Mat& v: T) { Mat w = mul(r.a, v); totalCandidates++; if (rightTerminal(w)) Tn.push_back(w); }
        sort(Tn.begin(),Tn.end()); size_t before=Tn.size(); Tn.erase(unique(Tn.begin(),Tn.end()),Tn.end());
        if (before!=Tn.size()) { fprintf(stderr,"duplicate candidates (coset decomposition not unique?)\n"); exit(1); }
        string Js; for(int j:Jp) Js.push_back('0'+j);
        printf("step %zu: J'={%s} cosets=%zu right-terminals(prev)=%zu candidates=%zu -> right-terminals=%zu\n", step, Js.c_str(), reps.size(), T.size(), reps.size()*T.size(), Tn.size());
        T.swap(Tn); J=Jp;
    }
    printf("total candidates tested: %lld\n", totalCandidates);
    // two-sided
    long long two=0, comm=0, noncomm=0; vector<string> words;
    for (const Mat& m: T) { Mat inv=inverse(m); if (rightTerminal(inv)) { two++; if (isCommutingProduct(m)) comm++; else { noncomm++; words.push_back(word(m)+"  len="+to_string(length(m))+" "+descents(m,inv)); } } }
    printf("E%d two-sided terminals=%lld (commuting %lld, noncommuting %lld)\n", N, two, comm, noncomm);
    sort(words.begin(),words.end(),[](const string&a,const string&b){return a.size()<b.size()||(a.size()==b.size()&&a<b);});
    for (auto& w: words) printf("  noncommuting terminal: %s\n", w.c_str());
}


// Bruhat test: x <= w iff (w==e ? x==e : with s in L(w): min(x,sx) <= sw)
static bool bruhatLeq(Mat x, Mat w) {
    while (true) {
        int s=-1; for (int i=0;i<N;i++) { Mat wi = inverse(w); if (isRightDescent(wi,i)) { s=i; break; } }
        if (s<0) return x==identity();
        Mat xi = inverse(x);
        if (isRightDescent(xi,s)) leftMul(x,s); // x -> sx when sx<x
        leftMul(w,s);
    }
}
static void eligible(const vector<string>& words) {
    // FC catalogue
    vector<Mat> all; unordered_set<Mat,MatHash> prev; prev.insert(identity()); all.push_back(identity());
    while (!prev.empty()) {
        unordered_set<Mat,MatHash> cand;
        for (const Mat& m: prev) for (int i=0;i<N;i++) if (!isRightDescent(m,i)) { Mat t=m; rightMul(t,i); cand.insert(t); }
        unordered_set<Mat,MatHash> nxt;
        for (const Mat& u: cand) {
            vector<int> R; for(int s=0;s<N;s++) if(isRightDescent(u,s)) R.push_back(s);
            bool ok=true;
            for (size_t i=0;i<R.size()&&ok;i++) for (size_t j=i+1;j<R.size();j++) if (ADJ[R[i]][R[j]]) {ok=false;break;}
            for (size_t i=0;i<R.size()&&ok;i++) { Mat t=u; rightMul(t,R[i]); if (!prev.count(t)) ok=false; }
            if (ok) nxt.insert(u);
        }
        for (auto& u: nxt) all.push_back(u); prev.swap(nxt);
    }
    printf("FC catalogue size %zu\n", all.size());
    for (const string& w: words) {
        Mat b=identity(); for (char ch: w) rightMul(b, ch-'0'); Mat bi=inverse(b);
        int cnt=0; string found;
        for (const Mat& x: all) {
            Mat xi=inverse(x); bool ok=true;
            for (int s=0;s<N&&ok;s++) { if (isRightDescent(b,s) && !isRightDescent(x,s)) ok=false; if (isRightDescent(bi,s) && !isRightDescent(xi,s)) ok=false; }
            if (!ok) continue;
            if (!bruhatLeq(x,b)) continue;
            cnt++; found += " " + word(x) + "(len " + to_string(length(x)) + ")";
        }
        printf("terminal %s (len %d): eligible FC lower endpoints: %d ->%s\n", w.c_str(), length(b), cnt, found.c_str());
    }
}


static void starChecks() {
    vector<Mat> all; unordered_set<Mat,MatHash> prev; prev.insert(identity()); all.push_back(identity());
    while (!prev.empty()) {
        unordered_set<Mat,MatHash> cand;
        for (const Mat& m: prev) for (int i=0;i<N;i++) if (!isRightDescent(m,i)) { Mat t=m; rightMul(t,i); cand.insert(t); }
        unordered_set<Mat,MatHash> nxt;
        for (const Mat& u: cand) {
            vector<int> R; for(int s=0;s<N;s++) if(isRightDescent(u,s)) R.push_back(s);
            bool ok=true;
            for (size_t i=0;i<R.size()&&ok;i++) for (size_t j=i+1;j<R.size();j++) if (ADJ[R[i]][R[j]]) {ok=false;break;}
            for (size_t i=0;i<R.size()&&ok;i++) { Mat t=u; rightMul(t,R[i]); if (!prev.count(t)) ok=false; }
            if (ok) nxt.insert(u);
        }
        for (auto& u: nxt) all.push_back(u); prev.swap(nxt);
    }
    unordered_set<Mat,MatHash> fcset(all.begin(), all.end());
    long long checks=0, failures=0;
    for (const Mat& x: all) {
        for (int side=0; side<2; side++) {
            Mat y = side==0 ? x : inverse(x); // right stars of x = right stars of x^{-1} inverted
            for (int s=0;s<N;s++) for (int t=s+1;t<N;t++) if (ADJ[s][t]) {
                int ds = isRightDescent(y,s), dt = isRightDescent(y,t);
                if (ds+dt!=1) continue;
                int found=0; Mat img;
                for (int u: {s,t}) { Mat z=y; rightMul(z,u); int a=isRightDescent(z,s)+isRightDescent(z,t); if (a==1) { found++; img=z; } }
                if (found!=1) { failures++; continue; }
                Mat imgx = side==0 ? img : inverse(img);
                if (!fcset.count(imgx)) failures++;
                checks++;
            }
        }
    }
    printf("E%d: FC=%zu star checks=%lld failures=%lld\n", N, all.size(), checks, failures);
}

// ---------- check a given word ----------
static void checkWord(const string& w) {
    Mat m=identity(); int len=0;
    for (char ch: w) { int i=ch-'0'; if (isRightDescent(m,i)) { printf("word %s NOT reduced at position %d\n", w.c_str(), len); } rightMul(m,i); len++; }
    Mat inv=inverse(m);
    printf("word %s: canonical=%s length=%d (by stripping: %d) supp=%s %s rightTerminal=%d leftTerminal=%d commutingProduct=%d\n", w.c_str(), word(m).c_str(), len, length(m), support(w).c_str(), descents(m,inv).c_str(), (int)rightTerminal(m), (int)rightTerminal(inv), (int)isCommutingProduct(m));
}

int main(int argc, char** argv) {
    if (argc<3) { fprintf(stderr,"usage: en <n> full|fc|prune <order>|word <w>...\n"); return 1; }
    int n=atoi(argv[1]); setup(n); setupInverse();
    printf("det(C)=%lld\n", detC);
    string mode=argv[2];
    if (mode=="full") fullEnumerate(true);
    else if (mode=="fc") fcEnumerate();
    else if (mode=="prune") { vector<int> order; for (char ch: string(argv[3])) order.push_back(ch-'0'); pruning(order); }
    else if (mode=="eligible") { vector<string> ws; for (int i=3;i<argc;i++) ws.push_back(argv[i]); eligible(ws); }
    else if (mode=="stars") starChecks();
    else if (mode=="word") { for (int i=3;i<argc;i++) checkWord(argv[i]); }
    return 0;
}
