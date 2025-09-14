#pragma once

#include <iostream>
#include <type_traits>
#include <vector>
#include <string>

/// @brief 2次元ボードを扱うクラス
/// @details 範囲外アクセスを許容して操作を単純にする
/// @details 範囲外をgetすると不正値を返す
/// @details 範囲外のsetは何もしない
template<typename T, T invalid_value = -1>
class Board{
private:
    std::vector<std::vector<T>> board;
    size_t height,width;

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