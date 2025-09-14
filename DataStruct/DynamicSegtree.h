#pragma once
#include <vector>
#include<functional>

/// @brief 動的セグ木、インデックスが連続でなくても区間取得が可能
/// @tparam T 
template<typename T>
class dynamic_segtree{
  private:
    using F = std::function<T(T,T)>;
    F f;
    F g;
    T ei;
    struct node_t{
        int L,R;
        T val;
        node_t *ch[2];
        node_t(int L,int R,T val):L(L),R(R),val(val){
            ch[0]=NULL;
            ch[1]=NULL;
        }
    };
    node_t *root;
    node_t *update(node_t* nd,int k,T x){
        if(nd->L==nd->R){
            nd->val = g(nd->val,x);
            return nd;
        }
        // Rが奇数のとき、Lの方向に丸める
        int r_true = nd->R + nd->L;
        if(r_true < 0) r_true -= 1;
        r_true /= 2;
        if(k<=r_true){
            if(!nd->ch[0]) nd->ch[0] = new node_t(nd->L,r_true,ei);
            nd->ch[0] = update(nd->ch[0],k,x);
        }
        else{
            if(!nd->ch[1]) nd->ch[1] = new node_t(r_true+1,nd->R,ei);
            nd->ch[1] = update(nd->ch[1],k,x);
        }
        nd->val = f(nd->ch[0]?nd->ch[0]->val:ei,nd->ch[1]?nd->ch[1]->val:ei);
        return nd;
    }
    T query(node_t *nd,int l,int r){
        if(nd->R<l||r<nd->L) return ei;
        if(l<=nd->L&&nd->R<=r) return nd->val;
        return f(nd->ch[0]?query(nd->ch[0],l,r):ei,nd->ch[1]?query(nd->ch[1],l,r):ei);
    }
  public:
    
    dynamic_segtree(F f,T ei, int range_l, int range_r):f(f),ei(ei){
        root = new node_t(range_l, range_r, ei);
    }
    ~dynamic_segtree(){delete root;root = nullptr;}
    /// @brief 一点更新
    /// @param k : index
    /// @param x  : 作用モノイド
    /// @param fc : 作用関数 
    void update(int k,T x,F fc){
        g=fc;
        update(root,k,x);
    }

    /// @brief 区間取得(閉区間)
    /// @param l : 左端
    /// @param r : 右端
    /// @return l以上r以下のインデックスにおける累積
    T query(int l,int r){
        return query(root,l,r);
    }
};