#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <stdexcept>

namespace kt {

    template <typename T>
    class vector2D {
    private:
        size_t rows_ ;
        size_t cols_ ;
        std :: vector <T > data_ ;

    public:
        vector2D (size_t rows , size_t cols , T val =0)
            : rows_( rows ) , cols_( cols ) , data_( rows * cols , val ) {}

        // Access element using [] operator
        auto* operator[](size_t row) {
            return & data_[row * cols_];
        }

        // Return a pointer to the first element
        T* data() {
            return data_.data () ;
        }

        // Other methods similar to std::vector
        void push_back(const T& value) {
            data_.push_back (value);
            cols_ ++;
        }

        size_t size() const {
            return data_.size () ;
        }

        size_t rows() const {
            return rows_;
        }

        size_t cols() const {
            return cols_;
        }
    };

}

#endif // VECTOR_H
