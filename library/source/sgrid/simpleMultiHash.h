#ifndef SIMPLEMULTIHASH_H
#define SIMPLEMULTIHASH_H


#include <math.h>
#include <cstring>
using namespace std;


/**
 * @brief A simple vector implementation so as not to rely on the standard library implementation.
 * This implementation has some of the same methods as the standard library.
 * 
 * @tparam T 
 */
template<typename T>
class SimpleVector{
private:
    T* _data;
    /// @brief reserved capacity for this vector. Always bigger or euqal to size
    size_t _capacity;
    /// @brief  current number of elements in the Vector
    size_t _size;
public:
    SimpleVector(size_t initSize = 0):_capacity(initSize+10),_size(initSize){
        // cout << "constructor" <<endl;
        _data = new T[_capacity]();
        // _data = (T*) calloc(_capacity,sizeof(T));
    }

    /**
     * @brief Destroy the Simple Vector object
     */
    ~SimpleVector(){
        // cout << "destructor" <<endl;
        delete [] _data;
        // free(_data);
    }

    /**
     * @brief copy constructor.
     * 
     * @param other
     */
    SimpleVector(const SimpleVector &other){
        // cout << "copy constructor" <<endl;
        _capacity = other._capacity;
        _size=other._size;
        _data = new T[_capacity]();
        // _data = (T*) calloc(_capacity,sizeof(T));

        memcpy(_data,other._data,_size*sizeof(T));
    }

    /**
     * @brief move constructor
     * 
     * @param other 
     */
    SimpleVector(SimpleVector&& other) noexcept 
    //: dataTableVector(std::exchange(other.dataTableVector, nullptr)) {} //TODO cpp14
    {
        // cout << "move constructor" <<endl;
        _data=other._data;
        other._data=nullptr;
        _size = other._size;
        _capacity =other._capacity;
    }
    /**
     * @brief copy assignment
     * 
     * @param other 
     * @return SimpleVector& 
     */
    SimpleVector& operator=(const SimpleVector& other) 
    {
        cout << "copy assignment" <<endl;
        if (this == &other) return *this;
        return *this = SimpleVector(other);
        // _capacity = other._capacity;
        // _data = new T[_capacity];
        // _size = other._size;
        // return *this;
    }
    /**
     * @brief move assignment
     * 
     * @param other 
     * @return SimpleVector& 
     */
    SimpleVector& operator=(SimpleVector&& other) noexcept 
    {
        // cout << "move assignment" <<endl;
        std::swap(_data, other._data);
        _size = other._size;
        _capacity = other._capacity;
        return *this;
    }

    T& operator [] (size_t i){ return at(i); _data[i]; }
    T& operator [] (size_t i) const { return at(i); _data[i]; }
    /**
     * @brief Array access with error checking
     * @throw out_of_range
     * 
     * @param i 
     * @return T& 
     */
    T& at(size_t i) {if(i>=_size)throw out_of_range("Illegal vector access"); return _data[i];}
    T& at(size_t i) const {if(i>=_size) {
        cout << "Cap: " << _capacity << "Size: " << _size << ", Idx: " << i <<endl;
        throw out_of_range("Illegal vector access2");
        } return _data[i];}

    void reserve(size_t newCapacity){
        // realloc() // HACK it is not safe to use realloc with NEW
        if(newCapacity <= _capacity) return;
        T* oldData = _data;
        _data = new T[newCapacity]();
        
        // _data = (T*) calloc(newCapacity,sizeof(T));
        memcpy(_data,oldData,_size * sizeof(T));
        delete [] oldData;
        //free(oldData);
        _capacity = newCapacity;
    }

    void push_back(const T& value){
        if(_size >= _capacity){reserve(_capacity<<1);} //TODO maybe other reserve strategy
        _data[_size] = value;
        _size++;
    }

    size_t size() const { return _size;}
    size_t capacity() const { return _capacity;}
    T* data() {return _data;}
    T* data() const {return _data;}
    

    T* begin(){return _data;};
    T* begin() const {return _data;};
    T* end(){return _data+_size;};
    T* end() const {return _data+_size;};
};


/**
 * @brief A simple multihash implementation.
 * Is able to store multiple elements with the same hash.
 * 
 * @tparam T 
 */
template<typename T, typename D>
class SimpleMultiHash{
    private:
    D depth;
     unsigned long hash(const T& point,const unsigned long tableSize) const {


        unsigned long value = point.getIndex(0);
        for (int d = 1; d < DimensionSparseGrid; ++d) {
            value = value + point.getIndex(d) * PrimeNumbers::getPrimeForHash(d);
        }

        value = value + depth.at(0) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid);
        for (int d = 1; d < DimensionSparseGrid; d++)
            value = value + depth.at(d) * PrimeNumbers::getPrimeForHash(DimensionSparseGrid + d);

        return value % tableSize;

    }

    SimpleVector<unsigned long> firstTable; // 0 means not used. Contains the pos+1 in secondTable  
    SimpleVector<unsigned long> secondTable; // 0  means unused. 1 means end of entries for same hash. Contains pos+2 of next entry
    SimpleVector<T> items; 
    //unsigned long (*hash)(const T& point,const unsigned long tableSize); // hash function


    unsigned long getEndOfHashgroup(const unsigned long index) const{
        unsigned long tmp = index;
        while (secondTable[tmp]>1)
        {
            tmp = secondTable[tmp]-2;
        }
        return tmp;
    }

/**
 * @brief Will try to find the given item in the Hashgroup
 * 
 * @param index 
 * @param item 
 * @param outFoundPos contains the position at which the item was found. If the item was not found returns the last position of the Hashgroup
 * @return true 
 * @return false 
 */
    bool checkIfItemInHashgroup(const unsigned long index,const T& item, unsigned long& outFoundPos) const {
        size_t tmp = index;
        while (secondTable[tmp]>=1)
        {
            if( items[tmp] == item ){
                outFoundPos = tmp;
                return true;
            }
            if(secondTable[tmp]==1) break;
            tmp = secondTable[tmp]-2;
        }
        outFoundPos = tmp;
        return false;
    }

    void checkIncreaseSecondTableSize() {
        if(secondTable.capacity() == secondTable.size()){
            secondTable.reserve(secondTable.size()*2);
            items.reserve(items.size()*2);
        }
    }

    public:

    unsigned long size() const {
        return secondTable.size();
    }
    /**
     * @brief Construct a new Simple Multi Hash object
     * 
     * @param hash 
     * @param firstTableSize The size of the first Table
     * @param secondTableInitSize The initially reserved space of the second table. May expand if set to low but will not shrink.
     */
    SimpleMultiHash(D depth_ , const unsigned long firstTableSize = 100, const unsigned long secondTableInitSize=300):firstTable(firstTableSize)/*,secondTable(secondTableSize)*/{
        depth=depth_;
        secondTable.reserve(secondTableInitSize);
        items.reserve(secondTableInitSize);
    }

    /**
     * @brief Will add the given point to the MultiHash. If the same point already exists it will not be added
     * 
     * @param point 
     * @return true 
     * @return false 
     */
    bool addPoint(const T& point){
        unsigned long pos = hash(point,firstTable.size());
        if(firstTable[pos]==0){
            firstTable[pos]=secondTable.size()+1;
            // secondTable.push_back(1);
            // return true;
        }else{
            unsigned long lastPos;
            if(checkIfItemInHashgroup(firstTable[pos]-1,point,lastPos)){
                return false;
            }
            secondTable[lastPos] = secondTable.size()+2;
            checkIncreaseSecondTableSize();
        }
        secondTable.push_back(1);
        items.push_back(point);
        return true;
    }

    // bool addPoint(const T& point, unsigned long& pos){
    //     size_t pos = hash(point,firstTable.size());
    //     if(firstTable[pos]==0){
    //         firstTable[pos]=secondTable.size()+1;
    //         // secondTable.push_back(1);
    //         // return true;
    //     }else{
    //         size_t lastPos;
    //         if(checkIfItemInHashgroup(firstTable[pos]-1,point,&lastPos)){
    //             return false;
    //         }
    //         secondTable[lastPos] = secondTable.size()+2;
    //         checkIncreaseSecondTableSize();
    //     }
    //     secondTable.push_back(1);
    //     items.push_back(point);
    //     pos = items.size()-1;
    //     return true;
    // }

    // bool addPointAtPos(const T& point, unsigned long pos){
    //     size_t pos = hash(point,firstTable.size());
    //     if(firstTable[pos]==0){
    //         firstTable[pos]=secondTable.size()+1;
    //         // secondTable.push_back(1);
    //         // return true;
    //     }else{
    //         size_t lastPos;
    //         if(checkIfItemInHashgroup(firstTable[pos]-1,point,&lastPos)){
    //             return false;
    //         }
    //         secondTable[lastPos] = secondTable.size()+2;
    //         checkIncreaseSecondTableSize();
    //     }
    //     secondTable.push_back(1);
    //     items.push_back(point);
    //     return true;
    // }
    /**
     * @brief Will add the given point to the MultiHash. If the point already exists inside the Hashmap the behaviour is undefined!!!!
     * Do NOT use this method if you are not absolutly sure that the given point does not already exist. Use addPoint() instead.
     * @param point 
     */
    void addPoint_Unique(const T& point){
        unsigned long pos = hash(point,firstTable.size());
        if(firstTable[pos]==0){
            firstTable[pos]=secondTable.size()+1;
        }else{
            unsigned long lastPos = getEndOfHashgroup(firstTable[pos]-1);
            secondTable[lastPos] = secondTable.size()+2;
            checkIncreaseSecondTableSize();
        }
        secondTable.push_back(1);
        items.push_back(point);
    }

    bool occupied(unsigned long& foundPos, const T& pointToSearch) const {
        unsigned long pos = hash(pointToSearch,firstTable.size());
        if(firstTable[pos]==0)
            return false;
        if(checkIfItemInHashgroup(firstTable[pos]-1,pointToSearch,foundPos)){
            return true;
        }
        return false;

    }

    T getIndexOfTable(unsigned long i) const {
        return items[i];
    }

    void print() const{
        cout << "First Table: "<<endl;
        for (auto &&i : firstTable)
        {
            cout << i << ", " <<endl;
        }
        
        // cout << "Second Table: "<<endl;
        // for (auto &&i : secondTable)
        // {
        //     cout << i << ", " <<endl;
        // }
    }

    unsigned long getSizeOfHashgroup(const unsigned long index) const{
        unsigned long tmp = index;
        unsigned long size = 1;
        while (secondTable[tmp]>1)
        {
            size++;
            tmp = secondTable[tmp]-2;
        }
        return size;
    }

    long double getSTDOfHash(){
        unsigned long sum =0;
        for (auto &&i : firstTable)
        {   
            if(i==0) continue;
            unsigned long numOfItems = getSizeOfHashgroup(i-1);
            sum += numOfItems;
        }
        long double mean = sum/firstTable.size();
        long double std= 0;
        // calc std
        for (auto &&i : firstTable)
        {   
            if(i==0) continue;
            unsigned long numOfItems = getSizeOfHashgroup(i-1);
            std += (numOfItems-mean)*(numOfItems-mean);
        }
        std /= firstTable.size() -1;
        return sqrt(std);
    }

    long double getRedDragonHashQuality(){
        long double res=0;
        for (auto &&i : firstTable)
        {   
            if(i==0) continue;
            unsigned long numOfItems = getSizeOfHashgroup(i-1);
            res += (numOfItems*(numOfItems+1)/2.0)/((double)(secondTable.size()+2*firstTable.size()-1));
        }
        return res;
    }

    void reorder(){

    }
};

#endif // SIMPLEMULTIHASH_H