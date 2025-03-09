#include<iostream>
#include<functional>

template<class Monoid, class EffectMonoid, bool IsCommutative = true>
class BinaryLinkedList{
private:
    using EffectFunc = std::function<Monoid(Monoid, EffectMonoid)>;
    struct Node{
        // ノード情報
        Node *nxt[2];
        Node *child_front[2];
        Monoid acc_value[2];
        bool rev;
        size_t size;
        size_t child_count;
        // 遅延評価
        EffectMonoid lazy;
        bool lazy_rev;
        explicit Node(const Monoid &value):
            nxt{nullptr, nullptr},
            child_front{nullptr, nullptr},
            acc_value{Monoid::id(), Monoid::id()},
            rev(false),
            size(1), // 末尾要素がいくつあるか
            child_count(0),
            lazy(EffectMonoid::id())
        {

        }

        bool debug_print(){
            std::cerr << "*****node start*****:" << this << std::endl;
			std::cerr << "[size, acc, rev acc, lazy, rev, child count ]\r\n=[";
			std::cerr << size << "," << acc_value[rev] << ", ";
			std::cerr << acc_value[!rev] << ", " << lazy << ", ";
            std::cerr << rev << ", " << child_count << "]" << std::endl; 
            std::vector<Node *> check;
			Node *child = child_front[rev];
			while(child){
                check.push_back(child->value);
                child = child->nxt[rev];
			}
			child = child_front[!rev];
            bool out = false;
            size_t i = check.size();
			while(child){
                if(check[--i] != child){
                    out = true;
                }
				child = child->nxt[!rev];
			}
            std::cout << "*****node end*****" << std::endl;
            return out;
        }
    };
    using NodePtr = Node *;
    Node *front_node;

    void pushdown(Node *node){
        if(node->lazy != EffectMonoid::id() || node->lazy_rev){
            // 遅延を各要素に伝播する
            Node *child = node->child_front[0];
            while(child){
                if(node->lazy != EffectMonoid::id()){
                    child->acc_value[0] = effect_func(child->acc_value[0], node->lazy);
                    if constexpr(IsCommutative){
                        child->acc_value[1] = child->acc_value[0];
                    }
                    else{
                        child->acc_value[1] = effect_func(child->acc_value[1], node->lazy);
                    }
                    child->lazy = EffectMonoid::op(child->lazy, node->lazy);
                }
                child->rev = child->rev ^ node->lazy_rev;
                child->lazy_rev = child->lazy_rev ^ node->lazy_rev;
                child = child->nxt[0];
            }
            node->lazy = EffectMonoid::id();
            node->lazy_rev = false;
        }
    }

    /// ノード内の要素を実際に反転することでrev変数をtoggleできる
    NodePtr toggle_rev(NodePtr node){
        std::vector<Node *> vec;
        Node  *child = node->child_front[0];
        Node *nxt;
        // 要素自身と要素の先頭へのリンクの反転で、反転ができる
        // N -f> 1 <p-n> 2 <p-n> 3 <b- N
        // N -f> 1 <n-p> 2 <n-p> 3 <b- N : element reverse
        // N -b> 1 <n-p> 2 <n-p> 3 <f- N : front reverse
        while(child){
            nxt = child->nxt[0];
            std::swap(child->nxt[0], child->nxt[1]);
            child = nxt;
        }
        std::swap(node->child_front[0], node->child_front[1]);
        // revのtoggleにより累積値も入れ替わる
        std::swap(node->acc_value[0], node->acc_value[1]);
        node->rev = !node->rev;
        return node;
    }

    void rev_accumulate(Node *node){
        Node *child = node->child_front[!node->rev];
        node->acc_value[!node->rev] = Monoid::id();
        while(child){
            node->acc_value[!node->rev] = Monoid::op(node->acc_value[!node->rev], child->acc_value[!child->rev]);
            child = child->nxt[!node->rev];
        }
    }

    void node_refresh(Node *node){
        node->size = 0;
        node->child_count = 0;
        Node *child = node->child_front[node->rev];
        if(!child){
            // 子がいないときは自身が末尾なので、累積値を行進しない
            return;
        }
        // 遅延評価で子を行進してから累積
        pushdown(node);
        node->acc_value[node->rev] = Monoid::id();
        while(child){
            node->size += child->size;
            node->child_count++;
            node->acc_value[node->rev] = Monoid::op(node->acc_value[node->rev], child->acc_value[child->rev]);
            child = child->nxt[node->rev];
        }
        // 逆方向の累積も行う
        if constexpr(IsCommutative){
            node->acc_value[!node->rev] = node->acc_value[node->rev];
        }
        else{
            rev_accumulate();
        }
    }

    /// @brief ノードを 1 | 2に分割
    /// @param node 分割対象のノード
    /// @param target_nxt ノード2の先頭になるノード
    /// @param root 分割対象の親ノード
    /// @return ノード1
    Node *split(Node *node, Node *target_nxt, Node *root){
        bool root_rev = root->rev;
        Node *target_prv = target_nxt->nxt[!root_rev];
        Node *nxt_node = new Node(Monoid::id());
        // ノード情報のコピー
        nxt_node->rev = node->rev;
        nxt_node->lazy = node->lazy;
        nxt_node->lazy_rev = node->lazy_rev;
        // 分割するノードのリンク張替え
        Node *jump_node = node->nxt[root_rev];
        nxt_node->nxt[root_rev] = jump_node;
        nxt_node->nxt[!root_rev] = node;
        node->nxt[root_rev] = nxt_node;
        if(jump_node){
            jump_node->nxt[!root_rev] = nxt_node;
        }
        // 子へのリンクの張替え
        nxt_node->child_front[nxt_node->rev] = target_nxt;
        nxt_node->child_front[!nxt_node->rev] = node->child_front[!node->rev];
        if(!target_prv){
            // 先頭で分割する場合
            node->child_front[node->rev] = nullptr;
        }
        node->child_front[!node->rev] = target_prv;
        // 子同士のリンク張替え
        if(target_prv){
            target_prv->nxt[node->rev] = nullptr;
        }
        target_nxt->nxt[!nxt_node->rev] = nullptr;
        // 親情報の更新
        if(root->child_front[!root_rev] == node){
            // 末尾を分割した場合ノード2を末尾とする
            root->child_front[!root_rev] = nxt_node;
        }
        root->child_count++;

        node_refresh(node);
        node_refresh(nxt_node);
        return node;
    }

    /// @brief 対象のノードと次のノードを連結する
    ///         親を超えて連結することはない
    /// @param node 対象ノード
    /// @param root 対象ノードの親
    /// @return 連結後のノード
    Node *connect(Node *node, Node *root){
        bool root_rev = root->rev;
        Node *nxt_node = node->nxt[root_rev];
        if(!nxt_node){
            return node;
        }
        // 向きを揃える
        if(node->rev != nxt_node->rev){
            toggle_rev(node);
        }
        // 遅延は接続前に伝播しておく
        node = pushdown(node);
        nxt_node = pushdown(nxt_node);

        // 子同士のリンク張替え
        Node *target_prv = node->child_front[!node->rev];
        Node *target_nxt = nxt_node->child_front[nxt_node->rev];
        if(target_prv){
            target_prv->nxt[node->rev] = target_nxt;
        }
        if(target_nxt){
            target_nxt->nxt[!nxt_node->rev] = target_prv;
        }
        // ノードから子へのリンク張替え
        if(nxt_node->child_count){
            node->child_front[!node->rev] = nxt_node->child_front[!nxt_node->rev];
            if(!node->child_count){
                node->child_front[node->rev] = target_nxt;
            }
        }
        // ノードのリンク張りかえ
        Node *jump_node = nxt_node->nxt[root_rev];
        node->nxt[root_rev] = jump_node;
        if(jump_node){
            jump_node->nxt[!root_rev] = node;
        }
        
        // 親からのリンクの更新
        if(nxt_node == root->child_front[!root_rev]){
            root->child_front[!root_rev] = node;
        }

        // ノード情報の更新
        node->size += nxt_node->size;
        node->child_count += nxt_node->child_count;
        node->acc_value[root_rev] = Monoid::op(node->acc_value[root_rev], nxt_node->acc_value[root_rev]);
        node->acc_value[!root_rev] = Monoid::op(nxt_node->acc_value[!root_rev], node->acc_value[!root_rev]);
        delete nxt_node;

        return node;
    }

    Node *input(Node *root, size_t idx, const Monoid &value){
        if(root->size == 1 && root->child_count == 0){
            // 要素を表すノードのとき
            return root;
        }
        Node *child = root->child_front[root->rev];
        Node *target;
        while(child && idx <= child->size){
            target = child;
            idx -= child->size;
            child = child->nxt[root->rev];
        }
        bool last = false;
        if(!child){
            idx = target->size;
        }
        else{
            target = child;
        }
        child = input(target, idx, value);
        if(child){

        }
    }

public:
    explicit BinaryLinkedList(EffectFunc q):
        front_node(new Node(Monoid::id()))
    {
        // 先頭ノードは木の根と同期する
        front_node->size = 0;
    }
    BinaryLinkedList(const std::vector<Monoid> &vec, EffectFunc q):
        BinaryLinkedList(q)
    {

    }
};

