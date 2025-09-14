#pragma once
#include <stack>
#include <vector>

/// @brief Heavy Light Decomposition
/// @details 頂点、辺のパス、部分木クエリをO(logN)で答える
/// @tparam Graph 隣接リスト表現のグラフ
class HLD{
  private:
    using Graph=std::vector<std::vector<int>>;
    int root,N;
    std::vector<int> size,par,Restore,newid,front,next,Depth,leaf;
    Graph graph;
    int dfs(int u){
        int cnt=1,MAX=0;
        for(auto v:graph[u]){
            if(v==par[u]) continue;
            par[v]=u;
            cnt+=dfs(v);
            if(size[v]>MAX){
                MAX=size[v];
                next[u]=v;
            }
        }
        return size[u]=cnt;
    }
    void hld(int r){
        std::stack<int> que;
        int id=0;
        que.push(r);
        while(!que.empty()){
            int u=que.top();que.pop();
            newid[u]=id;
            Restore[id++]=u;
            front[u]=u;
            while(next[u]!=-1){
                for(auto v:graph[u]){
                    if(v==par[u]) continue;
                    if(v!=next[u]){
                        Depth[v]=Depth[u]+1;
                        que.push(v);
                    }
                    else{
                        newid[v]=id;
                        Restore[id++]=v;
                        front[v]=front[u];
                        Depth[v]=Depth[u];
                    }
                }
                u=next[u];
            }
        }
    }
  public:
    HLD(Graph &g,int root=0):graph(g),root(root),N(g.size()),
    size(N),par(N,-1),Restore(N),newid(N),front(N),next(N,-1),Depth(N),leaf(N,-1){
        dfs(root);
        hld(root);
    }
    int lca(int u,int v)const{
        while(front[u]!=front[v]){
            if(Depth[u]>Depth[v]) std::swap(u,v);
            v = par[front[v]];
        }
        return newid[u]<newid[v]?u:v;
    }
    /// @brief u,v間のパスにおける連続したノードid頂点組クエリ
    /// 計算量 : O(logV)
    /// @param u パス端点にある頂点u
    /// @param v パスの端点にある頂点v
    /// @return パス(u,v)上にある頂点が全て含まれる、区間[l,r]の集合
    std::vector<std::pair<int,int>> vquery(int u,int v)const{
        std::vector<std::pair<int,int>> res;
        int rt=lca(u,v);
        while(Depth[rt]<Depth[u]){
            res.emplace_back(newid[front[u]],newid[u]);
            u = par[front[u]];
        }
        while(Depth[rt]<Depth[v]){
            res.emplace_back(newid[front[v]],newid[v]);
            v = par[front[v]];
        }
        res.emplace_back(newid[rt],std::max(newid[u],newid[v]));
        return res;
    }
    /// @brief u,v間のパスにおける連続した辺id頂点組クエリ]
    /// 計算量O(logV)
    /// @param u パス端点にある頂点u
    /// @param v パスの端点にある頂点v
    /// @return パス(u,v)上にある辺のidが全て含まれる、区間[l,r]の集合
    std::vector<std::pair<int,int>> equery(int u,int v)const{
        //頂点idから親に向かう辺の番号をid-1とする
        std::vector<std::pair<int,int>> res;
        int rt=lca(u,v);
        while(Depth[rt]<Depth[u]){
            res.emplace_back(newid[front[u]]-1,newid[u]-1);
            u = par[front[u]];
        }
        while(Depth[rt]<Depth[v]){
            res.emplace_back(newid[front[v]]-1,newid[v]-1);
            v = par[front[v]];
        }
        int R = std::max(newid[u],newid[v]);
        if(newid[rt]!=R) res.emplace_back(newid[rt],R-1);
        return res;
    }
    /// @brief 頂点uの部分木情報
    /// 計算量償却O(1)
    /// @param u 頂点
    /// @return 頂点uを根とする部分木の区間[l, r]
    std::pair<int,int> tquery(int u){
        if(leaf[u]!=-1) return std::make_pair(newid[u],leaf[u]);
        leaf[u]=newid[u];
        for(auto v:graph[u]){
            if(v==par[u]) continue;
            tquery(v);
            leaf[u]=std::max(leaf[u],leaf[v]);
        }
        return std::make_pair(newid[u],leaf[u]);
    }
    /// @brief  頂点uのid、頂点uから根に向かう辺はid-1
    int id(int u)const{return newid[u];}
    // idが示す頂点の元の番号
    int restore(int ID)const{return Restore[ID];}
    /// @brief 頂点uの深さ
    int depth(int u)const{return Depth[u];}
    // 頂点uの親
    int parent(int u)const{return par[u];}
};