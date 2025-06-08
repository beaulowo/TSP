#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <chrono> // timing

using namespace std;
using Coord = pair<double,double>;
const int INF = numeric_limits<int>::max();
static const int SKIP_PERCENT = 10;  // skip 10% of nodes

vector<Coord> read_tsp(const string &f) {
    ifstream in(f);
    if(!in) { cerr<<"Error opening "<<f<<"\n"; exit(1); }
    vector<Coord> C; string line;
    bool sec=false;
    while(getline(in,line)){
        if(!sec){ if(line.find("NODE_COORD_SECTION")!=string::npos) sec=true; continue; }
        if(line=="EOF"||line=="END") break;
        istringstream ss(line);
        int idx; double x,y;
        if(ss>>idx>>x>>y) C.emplace_back(x,y);
    }
    return C;
}

vector<vector<int>> make_dist(const vector<Coord> &C) {
    int n=C.size();
    vector<vector<int>> D(n, vector<int>(n, INF));
    for(int i=0;i<n;++i){ D[i][i]=0;
        for(int j=i+1;j<n;++j){
            double dx=C[i].first-C[j].first;
            double dy=C[i].second-C[j].second;
            int d=int(round(hypot(dx,dy)));
            D[i][j]=D[j][i]=d;
        }
    }
    return D;
}

vector<int> prim_mst(const vector<vector<int>> &D) {
    int n=D.size();
    vector<int> key(n,INF), par(n,-1);
    vector<bool> inM(n,false);
    key[0]=0;
    priority_queue<pair<int,int>, vector<pair<int,int>>, greater<>> pq;
    pq.emplace(0,0);
    while(!pq.empty()){
        auto [w,u]=pq.top(); pq.pop();
        if(inM[u]) continue;
        inM[u]=true;
        for(int v=0;v<n;++v)
            if(!inM[v] && D[u][v]<key[v]){
                key[v]=D[u][v]; par[v]=u;
                pq.emplace(key[v],v);
            }
    }
    return par;
}

vector<vector<int>> make_adj(const vector<int> &par) {
    int n=par.size();
    vector<vector<int>> A(n);
    for(int v=1;v<n;++v){ int u=par[v]; A[u].push_back(v); A[v].push_back(u);}    
    return A;
}

void dfs(int u, const vector<vector<int>> &A, vector<bool> &vis, vector<int> &tour) {
    vis[u]=true; tour.push_back(u);
    for(int v:A[u]) if(!vis[v]) dfs(v,A,vis,tour);
}

vector<int> random_skip_reinsert(const vector<int> &init, const vector<vector<int>> &D) {
    int n=init.size();
    int k=max(1, n*SKIP_PERCENT/100);
    vector<int> nodes(init.begin()+1, init.end());
    mt19937 gen(random_device{}()); shuffle(nodes.begin(),nodes.end(),gen);
    unordered_set<int> skip(nodes.begin(), nodes.begin()+k);
    vector<int> tour={init[0]};
    for(int v:init) if(v!=init[0] && !skip.count(v)) tour.push_back(v);
    for(int v:vector<int>(nodes.begin(),nodes.begin()+k)){
        int m=tour.size(), best=0, bd=INF;
        for(int i=0;i<m;++i){ int u=tour[i], w=(i+1<m?tour[i+1]:tour[0]);
            int d=D[u][v]+D[v][w]-D[u][w]; if(d<bd){ bd=d; best=i; }}
        if(best+1<m) tour.insert(tour.begin()+best+1,v);
        else tour.push_back(v);
    }
    return tour;
}

int main(int argc, char**argv){
    if(argc<2){ cerr<<"Usage: "<<argv[0]<<" file.tsp\n"; return 1; }
    auto C=read_tsp(argv[1]); int n=C.size(); if(n<2){ cerr<<"No coords\n"; return 1; }
    auto D=make_dist(C);
    auto par=prim_mst(D); auto A=make_adj(par);
    vector<bool> vis(n,false); vector<int> init;
    dfs(0,A,vis,init); init.push_back(0);
    
    auto t0=chrono::high_resolution_clock::now();
    auto tour=random_skip_reinsert(init,D);
    auto t1=chrono::high_resolution_clock::now();
    double tm=chrono::duration<double>(t1-t0).count();

    long long cost=0;
    for(int i=0;i<n;++i){ int u=tour[i], v=(i+1<n?tour[i+1]:tour[0]); cost+=D[u][v]; }
    
    cout<<"RS-MST cost: "<<cost<<"   time: "<<tm<<"s\n";
    return 0;
}
