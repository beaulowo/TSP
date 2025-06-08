#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>               // for timing

using namespace std;
using Coord = pair<double,double>;
const int INF = numeric_limits<int>::max();

// read TSPLIB .tsp
vector<Coord> read_tsp(const string &f){
    ifstream in(f);
    if(!in){ cerr<<"Error opening "<<f<<"\n"; exit(1); }
    vector<Coord> C; string line;
    bool sec=false;
    while(getline(in,line)){
        if(!sec){
            if(line.find("NODE_COORD_SECTION")!=string::npos) sec=true;
            continue;
        }
        if(line=="EOF"||line=="END") break;
        istringstream ss(line);
        int i; double x,y;
        if(ss>>i>>x>>y) C.emplace_back(x,y);
    }
    return C;
}

// build integer Euclid matrix
vector<vector<int>> make_dist(const vector<Coord> &C){
    int n=C.size();
    vector<vector<int>> D(n,vector<int>(n,INF));
    for(int i=0;i<n;++i){
        D[i][i]=0;
        for(int j=i+1;j<n;++j){
            double dx=C[i].first-C[j].first;
            double dy=C[i].second-C[j].second;
            int d=int(round(hypot(dx,dy)));
            D[i][j]=D[j][i]=d;
        }
    }
    return D;
}

// Prim → parent[]
vector<int> prim(const vector<vector<int>> &D){
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
        for(int v=0;v<n;++v){
            if(!inM[v] && D[u][v]<key[v]){
                key[v]=D[u][v];
                par[v]=u;
                pq.emplace(key[v],v);
            }
        }
    }
    return par;
}

// parent[] → adjlist
vector<vector<int>> make_adj(const vector<int> &par){
    int n=par.size();
    vector<vector<int>> A(n);
    for(int v=1;v<n;++v){
        int u=par[v];
        A[u].push_back(v);
        A[v].push_back(u);
    }
    return A;
}

// preorder DFS
void dfs(int u, const vector<vector<int>> &A,
         vector<bool> &vis, vector<int> &tour){
    vis[u]=true;
    tour.push_back(u);
    for(int v:A[u]) if(!vis[v]) dfs(v,A,vis,tour);
}

int main(int argc, char**argv){
    if(argc<2){ cerr<<"Usage: "<<argv[0]<<" file.tsp\n"; return 1; }
    auto C = read_tsp(argv[1]);
    int n=C.size();
    if(n<2){ cerr<<"No coords\n"; return 1; }
    auto D = make_dist(C);

    // start timing
    auto t0 = chrono::high_resolution_clock::now();

    auto P  = prim(D);
    auto A  = make_adj(P);
    vector<bool> vis(n,false);
    vector<int> tour;
    dfs(0,A,vis,tour);
    tour.push_back(0);

    long long cost=0;
    for(int i=0;i<n;++i){
        int u=tour[i], v=tour[(i+1)%n];
        cost += D[u][v];
    }

    // end timing
    auto t1 = chrono::high_resolution_clock::now();
    double secs = chrono::duration<double>(t1-t0).count();

    cout<<"2-Approx (CLRS) cost: "<<cost
        <<"   time: "<<secs<<"s\n";
    return 0;
}
