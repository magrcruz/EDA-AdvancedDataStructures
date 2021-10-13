// Copyright

#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <stack>
#include <algorithm>
#include "Point.hpp"

using namespace std;

template <size_t N, typename ElemType>
struct KDTreeNode {
    typedef std::pair<Point<N>, ElemType> value_type;

    value_type value;
    KDTreeNode* ptrs[2];

    KDTreeNode(const value_type& value) {
        this->value = value;
        ptrs[0] = 0;
        ptrs[1] = 0;
    }
    KDTreeNode(const KDTreeNode& kdNode) {
        this->value = kdNode.value;
        ptrs[0] = 0;
        ptrs[1] = 0;
    }
};

// STACK
template <size_t N, typename ElemType>
struct StackNode {
    KDTreeNode<N, ElemType>* node = 0;
    int state = 0;
    StackNode(KDTreeNode<N, ElemType>* n) :node(n) {};
};
// END-STACK

// QUEUE
template <size_t N, typename ElemType>
struct CQueueNode{
    const KDTreeNode<N, ElemType>* node = 0;
    double distancia;
    CQueueNode(const KDTreeNode<N, ElemType>* no, double dis) {
        node = no;
        distancia = dis;
    }
};

template <size_t N, typename ElemType>
struct CQueue {
    typedef CQueueNode<N, ElemType> CQNode;
    int limite;
    vector<CQNode> vec;
    CQueue(int lim) : limite(lim) {};
    
    void insertar(const KDTreeNode<N, ElemType>* node, double dis) {
        CQNode aux(node, dis);
        if (vec.empty()) {
            vec.push_back(aux);
            return;
        } 
        if (vec.size() < limite) vec.push_back(aux);
        int i = vec.size() - 1;
        for (; i > 0 && dis < vec[i-1].distancia; i--)
            vec[i] = vec[i - 1];
        if(dis < vec[i].distancia) vec[i] = aux;
    }

};

//KD-TREE
template <size_t N, typename ElemType>
class KDTree {
private:
    size_t dimension_;
    size_t size_;

public:
    KDTreeNode<N, ElemType>* root;
    typedef std::pair<Point<N>, ElemType> value_type;

    KDTree();

    ~KDTree();

    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);

    size_t dimension() const;

    size_t size() const;
    bool empty() const;

    bool contains(const Point<N>& pt) const;

    void insert(const Point<N>& pt, const ElemType& value);
    void insert(const value_type& value);

    ElemType& operator[](const Point<N>& pt);

    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;

    ElemType knn_value(const Point<N>& key, size_t k) const;

    std::vector<ElemType> knn_query(const Point<N>& key, size_t k) const;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {//
    dimension_ = N;
    root = 0;
    size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {//
    stack<StackNode<N, ElemType>> st;
    if (root) st.push(root);
    while (!st.empty()) {
        if (st.top().state == 0) {
            if (st.top().node) {
                st.top().state++;
                st.push(st.top().node->ptrs[0]);
            }
            else st.pop();
        }
        else if (st.top().state == 1) {
            st.top().state++;
            st.push(st.top().node->ptrs[1]);
        }
        else {
            delete st.top().node;
            st.pop();
        }
    }
    root = 0;
}

template <size_t N, typename ElemType>
void copyR(KDTreeNode<N, ElemType>* o, KDTreeNode<N, ElemType>** d) {
    if (o) {
        *d = new KDTreeNode<N, ElemType>(*o);
        copyR(o->ptrs[0], &((*d)->ptrs[0]));
        copyR(o->ptrs[1], &((*d)->ptrs[1]));
    }
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {//
    this->dimension_ = rhs.dimension_;
    this->size_ = rhs.size_;
    copyR(rhs.root, &root);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {//
    this->dimension_ = rhs.dimension_;
    this->size_ = rhs.size_;
    copyR(rhs.root, &root);
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {//
    return dimension_;
}


template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {//
    if (!root) return true;
    return false;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {//
    KDTreeNode<N, ElemType>* ptr = root, * prev = 0;
    int d = 0;
    bool dir;
    while (ptr) {
        if (pt == ptr->value.first) return 1;
        if (pt[d] <= ptr->value.first[d]) dir = 0;//izq
        else dir = 1;//der
        prev = ptr;
        ptr = ptr->ptrs[dir];
        d++;
        if (d >= dimension_) d -= dimension_;
    }
    return false;
}


template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {//
    if (empty()) root = new KDTreeNode<N, ElemType>(make_pair(pt, value));
    else {
        KDTreeNode<N, ElemType>* ptr = root;
        KDTreeNode<N, ElemType>* prev = root;
        size_t d = 0;
        bool dir;
        while (ptr) {
            if (pt == ptr->value.first) {
                ptr->value.second = value;
                return;
            }
            if (pt[d] <= ptr->value.first[d]) dir = 0;//izq
            else dir = 1;//der
            prev = ptr;
            ptr = ptr->ptrs[dir];
            d++;
            if (d >= dimension_) d -= dimension_;
        }
        prev->ptrs[dir] = new KDTreeNode<N, ElemType>(make_pair(pt, value));
    }
    size_++;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {//
    if (!root) {
        insert(pt, ElemType());
        return static_cast<ElemType&>(root->value.second);
    }
    KDTreeNode<N, ElemType>* ptr = root, * prev = 0;
    int d = 0;
    bool dir;
    while (ptr) {
        if (pt == ptr->value.first)
            return static_cast<ElemType&>(ptr->value.second);
        if (pt[d] <= ptr->value.first[d]) dir = 0;//izq
        else dir = 1;//der
        prev = ptr;
        ptr = ptr->ptrs[dir];
        d++;
        if (d >= dimension_) d -= dimension_;
    }
    prev->ptrs[dir] = new KDTreeNode<N, ElemType>(make_pair(pt, ElemType()));
    size_++;
    return static_cast<ElemType&>(prev->ptrs[dir]->value.second);
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {//
    KDTreeNode<N, ElemType>* ptr = root;
    KDTreeNode<N, ElemType>* prev = root;
    size_t d = 0;
    bool dir;
    while (ptr) {
        if (pt == ptr->value.first)
            return static_cast<ElemType&>(ptr->value.second);
        if (pt[d] <= ptr->value.first[d]) dir = 0;//izq
        else dir = 1;//der
        prev = ptr;
        ptr = ptr->ptrs[dir];
        d++;
        if (d >= dimension_) d -= dimension_;
    }
    throw out_of_range("out_of_range");
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    KDTreeNode<N, ElemType>* ptr = root;
    KDTreeNode<N, ElemType>* prev = root;
    size_t d = 0;
    bool dir;
    while (ptr) {
        if (pt == ptr->value.first)
            return const_cast<ElemType&>(ptr->value.second);
        if (pt[d] <= ptr->value.first[d]) dir = 0;//izq
        else dir = 1;//der
        prev = ptr;
        ptr = ptr->ptrs[dir];
        d++;
        if (d >= dimension_) d -= dimension_;
    }
    throw out_of_range("out_of_range");
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N>& key, size_t k) const {//
    vector<ElemType> vec = knn_query(key, k);
    ElemType eMax, eActual;//elmento maximo y actual
    int rMax,rActual;//conteo de repetidos
    sort(vec.begin(), vec.end());

    rMax = rActual = 1;
    eMax = eActual = vec[0];
    for (int i = 1; i < vec.size(); i++, rActual++) {
        if (eActual != vec[i]) {
            rActual = 0;
            eActual = vec[i];
        }
        else if (rActual + 1 > rMax) {
            rMax = rActual + 1;
            eMax = eActual;
        }
    }
    return eMax;
}

template <size_t N, typename ElemType>
void knn_recursion(const Point<N>& pBuscado, const KDTreeNode<N, ElemType>* nActual, int dim, double &minDistance,CQueue<N,ElemType> &cola) {
    if (!nActual) return;
    dim %= N;
    double distancia = distance(pBuscado, nActual->value.first);
    cola.insertar(nActual, distancia);
    if (distancia <= minDistance) minDistance = distancia;
    double distancias[2] = { DBL_MAX, DBL_MAX };
    if (nActual->ptrs[0]) distancias[0] = distance(pBuscado, nActual->ptrs[0]->value.first);
    if (nActual->ptrs[1]) distancias[1] = distance(pBuscado, nActual->ptrs[1]->value.first);
    bool min = distancias[0] >= distancias[1];
    knn_recursion(pBuscado, nActual->ptrs[min], dim+1, minDistance, cola);
    if (nActual->value.first[dim] - pBuscado[dim] < distancia)//Cruza el cuadrante
        knn_recursion(pBuscado, nActual->ptrs[!min], dim+1, minDistance, cola);
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key,
    size_t k) const {//
    std::vector<ElemType> values;
    double minDis = DBL_MAX;
    CQueue<N, ElemType> cola(k);
    knn_recursion(key, root, 0, minDis,cola);
    for (int i = 0; i < cola.vec.size(); i++) values.push_back(cola.vec[i].node->value.second);
    return values;
}

#endif  // SRC_KDTREE_HPP_
