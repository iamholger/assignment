#include <vector>
#include <array>
#include <cassert>
#include <iostream>
#include <stdexcept>


class NotImplemented : public std::logic_error
{
    public:
        NotImplemented() : std::logic_error("Function not implemented") { };
};

template<size_t num_rows, size_t num_columns, typename T>
class Matrix;

template<size_t num_rows, typename T>
class Vector : public Matrix<num_rows, 1, T>
{
    public:
        using Matrix<num_rows, 1, T>::operator =;
};

template<size_t num_rows, size_t num_columns, typename T>
class Matrix
{
    public:

        Matrix() : nrows_(num_rows), ncols_(num_columns), current_(0) {}

        // A posteriori check to prevent usage of items with not enough inputs
        // provided on comma initialisation
        bool isCorrectlyInitialised() const {return current_ == this->size();}

        size_t nrows() const {return nrows_;}
        size_t ncols() const {return ncols_;}
        size_t size()  const {return ncols_*nrows_;}

        // Set element
        void setElement(size_t row, size_t col, T val)
        {
            if ( (row >= nrows_) or (col>=ncols_) ) throw std::out_of_range("Trying to set element out of range");
            data_[row][col] = val;
        }
        
        // Add value to element
        void addToElement(size_t row, size_t col, T val)
        {
            if ( (row >= nrows_) or (col>=ncols_) ) throw std::out_of_range("Trying to set element out of range");
            data_[row][col] += val;
        }
        
        // Get element
        T getElement(size_t row, size_t col) const
        {
            if ( (row >= nrows_) or (col>=ncols_) ) throw std::out_of_range("Trying to set element out of range");
            return data_[row][col];
        }

        // Determinant, closed expressions for 2D and 3D
        T determinant() const
        {
            if (ncols_ != nrows_) throw std::logic_error("Trying to compute determinant of non-square matrix");
            if ( nrows_>3 or nrows_<2) throw NotImplemented();
            T det(0);

            if (nrows_ == 2)
            {
                det = data_[0][0]*data_[1][1] - data_[0][1]*data_[1][0];
            }
            else if (nrows_ == 3)
            {
                det = data_[0][0]*data_[1][1]*data_[2][2] 
                    + data_[0][1]*data_[1][2]*data_[2][0]
                    + data_[0][2]*data_[1][0]*data_[2][1]
                    - data_[2][0]*data_[1][1]*data_[0][2]
                    - data_[2][1]*data_[1][2]*data_[0][0]
                    - data_[2][2]*data_[1][0]*data_[0][1];
            }
            else  std::logic_error("This should never be reached");

            return det;
        }

        // Co-factor computation of single matrix element Mij:
        //
        // cf_ij = (-1)^(i+j) * det(Minor_ij) with Minor_ij being Mij after removal
        // of row and and column j
        T cofactor(size_t row, size_t col) const
        {
            if ( (row >= nrows_) or (col>=ncols_) ) throw std::out_of_range("Trying to access element out of range");
            if (nrows_ !=ncols_)  throw std::logic_error("Cannot compute cofactor for non-square matrix");
            if (nrows_ != 3) throw NotImplemented();

            //  (-1)^(i+j)
            T cf = -1;
            if ( (row+col)%2 == 0 ) cf = 1;
            
            // Construct Minor
            Matrix<2,2,T> minor;
            minor.set_entries(0);

            int counter(0);
            for (int r=0;r<3;r++)
            for (int c=0;c<3;c++)
            {
                if (r!=row and c !=col)
                {
                    const size_t target_row = counter / 2;
                    const size_t target_col = counter % 2;
                    minor.setElement(target_row, target_col, data_[r][c]);
                    counter++;
                }
            }

            cf *= minor.determinant();

            return cf;
        }

        // Set all entries to val
        void set_entries(T val)
        {
            std::fill(*data_, *data_+this->size(),val);
            current_ = ncols_*nrows_;
        }

        // Return a new matrix of the same shape with all elements set to 0
        Matrix zeros_like() const
        {
            Matrix<num_rows,num_columns,T> out;
            out.set_entries(0);
            return out;
        } 


        // For priting to a stream
        friend std::ostream& operator<<(std::ostream& out, const Matrix& m)
        {
            if (not m.isCorrectlyInitialised()) throw std::logic_error("Matrix not properly initialised");
            out << "Matrix of shape " << m.nrows() << "x" << m.ncols() << "\n";
            for (int r=0;r<m.nrows();r++)
            {
                for (int c=0;c<m.ncols();c++)
                {
                    out << m.data_[r][c] << " ";
                }
                out << "\n";
            }
            out << "\n";
            return out;
        }
        
        
        // Helper struct for comma initialisation
        struct ReadComma
        {
            Matrix& m_; // The matrix we want to populate with data
            ReadComma(Matrix<num_rows,num_columns, T> & m) : m_(m) {}

            // overload operator, --- chaining operators 
            ReadComma operator,(const T x)
            {
                // Throw if too many elements provided by user
                if (m_.current_ >= m_.size()) throw std::logic_error("Trying to set more elements than possible");

                // Logic to figure out position in data array
                const size_t row = (m_.ncols()==1) ? m_.current_ : m_.current_ / m_.nrows();
                const size_t col = (m_.ncols()==1) ?            0: m_.current_ % m_.ncols();
#if DEBUG
                std::cout << "Input #" << m_.current_ << " sets element (" << row<<"," << col << ") to " << x <<"\n";
#endif
                m_.data_[row][col] = x;
                m_.current_++;
                return ReadComma(m_);
            }
        };
       
        // Overload operator= to allow for comma initialisation.
        ReadComma operator=(T x)
        {
            current_ = 0;    // Reset reader counter
            data_[0][0] = x; // Set first element directly
            current_++;
            return ReadComma(*this); // This will go and capture all the other inputs
        }
        
        // Assignment operator
        Matrix<num_rows, num_columns,T> & operator=(Matrix<num_rows, num_columns,T> other)
        {
            for (int r=0;r<other.nrows();r++)
            for (int c=0;c<other.ncols();c++)
            {
                this->setElement(r,c,other.getElement(r,c));
            }
            current_ = this->size(); // Set reader counter to signal correct initialisation
            return *this;
        }


        // Matrix vector multiplication Ci = Aij*Bj 
        friend Vector<num_columns,T> operator*(const Matrix<num_rows, num_columns,T> & mx, const Vector<num_rows, T> & vr)
        {
            Vector<num_columns,T> ret;
            ret.set_entries(0);
            
            for (int r=0;r<num_rows;   r++)
            for (int c=0;c<num_columns;c++) 
            {
                ret.addToElement(r, 0,  mx.getElement(r,c) * vr.getElement(c,0));
            }

            return ret;
        }

    private:
        T data_[num_rows][num_columns];
        const size_t nrows_;
        const size_t ncols_;
        size_t current_;

};

template<size_t num_rows, size_t num_columns, typename T>
Matrix<num_columns,num_rows, T> transpose(const Matrix<num_rows,num_columns, T> & mx)
{
    if (not mx.isCorrectlyInitialised()) throw std::logic_error("Matrix not properly initialised");
    Matrix<num_columns,num_rows, T> trp;
    trp.set_entries(0);

    for (int r=0;r<mx.nrows();r++)
    for (int c=0;c<mx.ncols();c++)
    {
        trp.setElement(c,r, mx.getElement(r,c));
    }

    return trp;
}

template<size_t num_rows, size_t num_columns, typename T>
Matrix<num_rows,num_columns, T> inverse(const Matrix<num_rows,num_columns, T> & mx)
{
    if (not mx.isCorrectlyInitialised()) throw std::logic_error("Matrix not properly initialised");
    if (mx.nrows() != mx.ncols())  throw std::domain_error("Matrix is not invertible");
    if (mx.nrows() != 3) throw NotImplemented();
    const T det = mx.determinant();
    if (det==0) throw std::domain_error("Matrix is not invertible");
    auto out = mx.zeros_like();


    Matrix<num_rows,num_columns, T> inv;
    inv.set_entries(1);

    for (int r=0;r<mx.nrows();r++)
    for (int c=0;c<mx.ncols();c++)  inv.setElement(r,c, mx.cofactor(r,c)/det);

    return transpose(inv);
}
