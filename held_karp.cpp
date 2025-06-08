#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>  // timing

using namespace std;
using Coord = pair<double,double>;
const int INF = numeric_limits<int>::max();

vector<Coord> read_tsp(const string &file) {
    ifstream in(file);
    if(!in) { cerr<<"Error opening "<<file<<"\n"; exit(1); }
    vector<Coord> coords;
    string line;
    bool section = false;
    while(getline(in,line)) {
        if(!section) {
            if(line.find("NODE_COORD_SECTION")!=string::npos) section = true;
            continue;
        }
        if(line=="EOF"||line=="END") break;
        istringstream ss(line);
        int idx; double x,y;
        if(ss>>idx>>x>>y) coords.emplace_back(x,y);
    }
    return coords;
}

vector<vector<int>> compute_dist(const vector<Coord> &C) {
    int n=C.size();
    vector<vector<int>> D(n, vector<int>(n, INF));
    for(int i=0;i<n;++i) {
        D[i][i]=0;
        for(int j=i+1;j<n;++j) {
            double dx=C[i].first-C[j].first;
            double dy=C[i].second-C[j].second;
            int d=int(round(hypot(dx,dy)));
            D[i][j]=D[j][i]=d;
        }
    }
    return D;
}

pair<int, vector<int>> held_karp(const vector<vector<int>> &dist) {
    int n = dist.size();
    int FULL = 1<<n;
    vector<vector<int>> dp(FULL, vector<int>(n, INF));
    vector<vector<int>> parent(FULL, vector<int>(n, -1));
    dp[1][0] = 0;
    for(int mask=0; mask<FULL; ++mask) {
        if(!(mask & 1)) continue;
        for(int u=0; u<n; ++u) {
            if(!(mask & (1<<u)) || dp[mask][u]==INF) continue;
            int rem = (FULL-1) ^ mask;
            for(int sub=rem; sub; sub &= sub-1) {
                int v = __builtin_ctz(sub);
                int nm = mask | (1<<v);
                int cost = dp[mask][u] + dist[u][v];
                if(cost < dp[nm][v]) {
                    dp[nm][v] = cost;
                    parent[nm][v] = u;
                }
            }
        }
    }
    int endMask = FULL-1, best=INF, last=0;
    for(int i=1;i<n;++i) {
        int c = dp[endMask][i] + dist[i][0];
        if(c<best) { best=c; last=i; }
    }
    vector<int> tour;
    int mask=endMask, cur=last;
    while(cur!=0) {
        tour.push_back(cur);
        int p = parent[mask][cur];
        mask ^= 1<<cur;
        cur = p;
    }
    tour.push_back(0);
    reverse(tour.begin(), tour.end());
    tour.push_back(0);
    return {best, tour};
}

int main(int argc, char **argv) {
    if(argc<2) { cerr<<"Usage: "<<argv[0]<<" file.tsp\n"; return 1; }
    auto coords = read_tsp(argv[1]);
    int n = coords.size();
    if(n<2) { cerr<<"No coordinates\n"; return 1; }
    auto dist = compute_dist(coords);

    auto t0 = chrono::high_resolution_clock::now();
    auto res = held_karp(dist);
    auto t1 = chrono::high_resolution_clock::now();
    double tm = chrono::duration<double>(t1-t0).count();

    int cost = res.first;
    cout<<"Held-Karp optimal cost: "<<cost
        <<"   time: "<<tm<<"s\n";
    return 0;
}
