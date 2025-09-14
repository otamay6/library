#pragma once

#include <iostream>
#include <type_traits>
#include <vector>
#include <string>
#include <functional>
#include <queue>

/// @brief 2次元ボードを扱うクラス
/// @details 範囲外アクセスを許容して操作を単純にする
/// @details 範囲外をgetすると不正値を返す
/// @details 範囲外のsetは何もしない
template<typename T, T invalid_value = -1>
class Board{
private:
    std::vector<std::vector<T>> board;
    size_t height,width;

    using SearchFunc = std::function<void(Board &, int, int)>;

    bool inside(int x, int y){
        return 0 <= x && x < (int)width && 0 <= y && y < (int)height;
    }
public:
    Board(size_t H, size_t W) : height(H), width(W){
        board.resize(H);
    }

    T get(int x, int y){
        if(inside(x, y)){
            return board[y][x];
        }
        return invalid_value;
    }

    void set(int x, int y, T value){
        if(inside(x,y)){
            board[y][x] = value;
        }
    }

    /// @brief 探索用の関数
    /// @param x,y : 探索初期位置
    /// @param search_list : [dx, dy]のリスト形式の移動候補
    /// @param func : 候補移動時に呼ばれる関数
    void search(int x, int y, const std::vector<std::array<int, 2>> &search_list, SearchFunc func){
        for(auto dir : search_list) {
            int dx = dir[0], dy = dir[1];
            func(*this, x + dx, y + dy);
        }
    }

    /// @brief 隣接する四方を探索するためのsearch関数ラッパ
    void search_neigher(int x, int y, SearchFunc func){
        search(x, y, {{-1,0}, {0, -1}, {1, 0}, {0, 1}}, func);
    }

    /// @brief 周囲8マスを探索するためのsearch関数ラッパ
    void search_round(int x, int y, SearchFunc func){
        search(x, y, 
            {
                {-1, -1}, {-1, 0}, {-1, 1}, 
                {0, -1}, {0, 1},
                {1, -1}, {1, 0}, {1, -1},
            }, func);
    }

    /// @brief BPFで探索する関数
    void search_bfs(int x, int y, SearchFunc func){
        std::queue<std::array<int , 2>> pos;
        std::vector<std::vector<bool>> searched(height, std::vector<bool>(width, false));
        pos.push({x, y});
        while(pos.size()){
            auto p = pos.front();
            pos.pop();
            int nx = p[0], ny = p[1];
            if(searched[ny][nx]) continue;
            func(board, x, y);
            search_neighber(x, y, [&pos](Board &, int x, int y){
                pos.push({x, y});
            });
            searched[ny][nx] = true;
        }
    };


    friend std::istream &operator >>(std::istream &is, Board &bd){
        for(size_t i = 0; i < bd.height; i++){
            std::vector<char> row;
            if constexpr (std::is_same_v<T, char>) {
                std::string s;
                is >> s;
                for(size_t j = 0; j < bd.width; j++){
                    row.push_back(s[j]);
                }
            }
            else {
                for(size_t j = 0; j < bd.width; j++){
                    T c;
                    is >> c;
                    row.push_back(c);
                }
            }
            bd.board[i] = row;
        }
        return is;
    }
};