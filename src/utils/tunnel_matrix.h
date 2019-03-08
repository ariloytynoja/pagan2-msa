
#ifndef TUNNEL_MATRIX_H
#define TUNNEL_MATRIX_H

#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <stdexcept>



//using namespace boost;


using namespace std;

/**
 * @brief Tunnel_slice is a column or a row of a Tunnel_matrix.
 *
 * Tunnel_slice is a container for a subset of the whole matrix's elements: a specific column or row.
 * Column slices contain the elements themselves and row slices only pointers the elements. This way
 * the elements are stored only once. Elements (or pointers) are stored in a one dimensional vector.
 * Element representing all elements outside the tunnel is saved as an pointer.
 *
 */
template <class Type> class Tunnel_slice
{
private:
    int total_length;
    int begin;
    int end;
    Type* outside_tunnel_entry;
    Type* entries;
    Type** entry_pointers;
    bool is_pointer_slice;

public:
    Tunnel_slice(int length, int tunnel_begin, int tunnel_end, Type* empty, bool use_pointers = false){
        is_pointer_slice = use_pointers;
        total_length = length;
        begin = tunnel_begin;
        end = tunnel_end;
        outside_tunnel_entry = empty;


        if(use_pointers){
            entry_pointers = new Type*[max(0, end - begin + 1)];
            entries = 0;
        }else{
            entries = new Type[max(0, end - begin + 1)];
            entry_pointers = 0;
        }

    }
    Tunnel_slice(int length, Type* empty, bool use_pointers = false){
        is_pointer_slice = use_pointers;
        total_length = length;
        begin = 0;
        end = -1;
        outside_tunnel_entry = empty;

    }
    Type& operator[](int index){
        return at(index);
    }
    Type& at(int index){
        if(index < begin || index > end){
            if(index < 0 || index > total_length){
                throw std::out_of_range("Tunnel_matrix: y-index out of bounds");
            }
            return *outside_tunnel_entry;

        }
        if(is_pointer_slice){
            return *entry_pointers[index - begin];
        }
        return entries[index - begin];

    }
    Type* ptrAt(int index){
        if(index < begin || index > end){
            if(index < 0 || index > total_length){
                throw std::out_of_range("Tunnel_matrix: pointer index out of bounds");
            }
            //TODDOreturn *outside_tunnel_entry;
            throw std::out_of_range("Tunnel_matrix: no pointers outside tunnel");
        }
        if(is_pointer_slice){
            return entry_pointers[index - begin];
        }
        return &entries[index - begin];
        //return entries.at(index - begin);
    }
    void setPtrAt(int index, Type* ptr){
        if(index < begin || index > end){
            if(index < 0 || index > total_length){
                throw std::out_of_range("Tunnel_matrix: pointer index out of bounds");
            }
            //TODDOreturn *outside_tunnel_entry;
            throw std::out_of_range("Tunnel_matrix: no pointers outside tunnel");
        }
        if(is_pointer_slice){
            entry_pointers[index - begin] = ptr;
        }else{
            throw std::logic_error("Tunnel_matrix: cannot set pointer of a non-pointer slice");
        }
    }
    int getBegin(){
        return begin;
    }
    int getEnd(){
        return end;
    }
    int getEntriesSize(){
        return end - begin + 1;
    }
    void setBegin(int index){
        begin = index;
    }
    void setEnd(int index){
        end = index;
    }
    int getTotalSize(){
        return total_length;
    }
    ~Tunnel_slice(){
        if(is_pointer_slice){
            delete []this->entry_pointers;
        }else{
            delete []this->entries;
        }
    }
};

/**
 * @brief Tunnel_matrix is a two dimensional data structure which consist of elements only inside a specific 'tunnel'.
 *
 * To reduce memory consumption only a specific tunnel inside NxM matrix is allocated. This is done by using N Tunnel_slices,
 * arrays that contains elements of a specific column (plus the starting and ending coordinates of the tunnel). M slices
 * of pointers are also included to allow accessing a specific row. One pointer is also reserved to represent all the elements
 * outside the tunnel.
 *
 * The tunnel is given as two vectors: upper and lower bounds (length N). The tunnel also has one restriction: it has to be
 * monotonically increasing. This allows a linear construction time.
 *
 */
template <class Type> class Tunnel_matrix
{
private:
    bool tunneled;
    int x_length;
    int y_length;
    vector<Tunnel_slice<Type>*> x_slices;
    vector<Tunnel_slice<Type>*> y_slices;
    Type* outside_tunnel_entry;
    vector<unsigned int> dims;

    void initMatrixSize(int dimx_length, int dimy_length){
        x_length = dimx_length;
        y_length = dimy_length;

        dims.push_back((unsigned int)x_length);
        dims.push_back((unsigned int)y_length);

    }
    void initTunnelEntries(vector<int>& x_starts, vector<int>& x_ends){

        tunneled = true;

        x_slices.reserve(x_length);
        y_slices.reserve(y_length);

        ///x_slices (with objects)
        for(int x = 0; x < x_length; x++){
            Tunnel_slice<Type>* x_slice = new Tunnel_slice<Type>(y_length, max(0, x_starts[x]), min(x_ends[x], y_length - 1), outside_tunnel_entry);
            x_slices.push_back(x_slice);
        }

        ///y_slices (with pointers to objects
        int x_start_index = 0;
        for(int y = 0; y < y_length; y++){


            int tunnel_width = 0;
            int tunnel_begin = 0;
            int tunnel_end = -1;
            for(int x = x_start_index; x < x_length;x++){
                Tunnel_slice<Type>& slx = slice_x(x);


                if(y >= slx.getBegin() && y <= slx.getEnd()){ //We are inside tunnel!

                    if(tunnel_width == 0){
                        x_start_index = x;
                        tunnel_begin = x;
                    }
                    tunnel_width++;
                    tunnel_end = x;
                }else if(y > slx.getEnd()){ //Haven't found tunnel yet..

                }else{  //Tunnel ended
                    break;
                }
            }

            Tunnel_slice<Type>* y_slice = new Tunnel_slice<Type>(x_length, tunnel_begin, tunnel_end, outside_tunnel_entry, true);

            for(int x = tunnel_begin; x <= tunnel_end; x++){
                Tunnel_slice<Type>& slx = slice_x(x);
                y_slice->setPtrAt(x, slx.ptrAt(y));
            }
            y_slices.push_back(y_slice);
        }
    }
    void initFullMatrixEntries(){

        x_slices.reserve(x_length);
        y_slices.reserve(y_length);

        ///x_slices (with objects)
        for(int x = 0; x < x_length; x++){
            Tunnel_slice<Type>* x_slice = new Tunnel_slice<Type>(y_length, 0, y_length - 1, outside_tunnel_entry);
            x_slices.push_back(x_slice);
        }

        ///y_slices (with pointers to objects)
        for(int y = 0; y < y_length; y++){
            Tunnel_slice<Type>* y_slice = new Tunnel_slice<Type>(x_length, 0, x_length - 1, outside_tunnel_entry, true);

            for(int x = 0; x < x_length;x++){
                Tunnel_slice<Type>& slx = slice_x(x);
                y_slice->setPtrAt(x, slx.ptrAt(y));
            }
            y_slices.push_back(y_slice);
        }

    }
    void initWithTunnelBounds(int dimx_length, int dimy_length, vector<int>& x_starts, vector<int>& x_ends){

        if(dimx_length < 0 || dimy_length < 0){
            throw std::length_error("Tunnel_matrix: invalid matrix size");
        }


        initMatrixSize(dimx_length, dimy_length);

        if(x_starts.size() == 0 && (int) x_ends.size() == 0){
            initFullMatrixEntries();
        }else if(dimx_length != (int) x_starts.size() || dimx_length != (int) x_ends.size()){
            throw std::length_error("Tunnel_matrix: invalid tunnel bounds lengths");
        }else{

            /*
            float avg_start = 0;
            float avg_end = 0;

            for(int i = 0; i < dimx_length; i++){
                avg_start += x_starts.at(i);
                avg_end += x_ends.at(i);
            }
            avg_start /= dimx_length;
            avg_end /= dimx_length;

            cout << "Avg start: " << avg_start << endl;
            cout << "Avg end: " << avg_end << endl;
            cout << "Tunnel size: " << 100* (avg_end - avg_start) / dimy_length << "%" << endl;
            */

            initTunnelEntries(x_starts, x_ends);
        }
    }
public:
    Tunnel_matrix(int dimx_length, int dimy_length) : outside_tunnel_entry(new Type){

        initMatrixSize(dimx_length, dimy_length);
        initFullMatrixEntries();
    }
    Tunnel_matrix(int dimx_length, int dimy_length, vector<int>& x_starts, vector<int>& x_ends) : outside_tunnel_entry(new Type){

        initWithTunnelBounds(dimx_length, dimy_length, x_starts, x_ends);

    }
    Tunnel_matrix(int dimx_length, int dimy_length, vector<int>& x_starts, vector<int>& x_ends, Type* empty){

        outside_tunnel_entry = empty;
        initWithTunnelBounds(dimx_length, dimy_length, x_starts, x_ends);
    }
    ~Tunnel_matrix(){
        for(int i = 0; i < (int) x_slices.size(); i++){
            delete x_slices.at(i);
            x_slices.at(i) = 0;
        }
        for(int i = 0; i < (int) y_slices.size(); i++){
            delete y_slices.at(i);
            y_slices.at(i) = 0;
        }
    }
    Type& at(int x, int y){
        return slice_x(x).at(y);
    }
    Tunnel_slice<Type>& operator[](int x){
        return slice_x(x);
    }
    Tunnel_slice<Type>& slice_x(int x){
        if(x < 0 || x >= x_length){
            throw std::out_of_range("Tunnel_matrix: x-index out of bounds");
        }
        return *x_slices.at(x);
    }
    Tunnel_slice<Type>& slice_y(int y){
        if(y < 0 || y >= y_length){
            throw std::out_of_range("Tunnel_matrix: y-index out of bounds");
        }
        return *y_slices.at(y);
    }
    int getSize_y(){
        return y_length;
    }
    int getSize_x(){
        return x_length;
    }
    vector<unsigned int>& shape(){
        return dims;
    }
};


#endif



