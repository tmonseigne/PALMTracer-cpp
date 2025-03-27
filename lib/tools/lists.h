#include "errc.h"
#include "generalParam.h"

#ifndef __PPLISTH__
#define __PPLISTH__

/*
 * Base class of Heap.
 */
template <class T>
class Theap
{
public :
	float key;
	T element;

	Theap() { }
	~Theap() { }
};

/**
 * The class <code>Heap</code> defines a heap of 
 * any primitive types (Uchar, Char, Slong, Float, ..).
 * In a heap, each element is indexed by a key value. The next element
 * to be popped is the element with the min key.
 * Element is inserted by <code>push()</code> and
 * removed by <code>pop()</code>.
 * To create the new Heap named "foo" with a maximum of 100 floats
 * use: <code>Heap<Float> foo(100)</code>.
 */
template <class T>
class Heap
{
private:
	int _size;
	int _nitem;
	Theap<T>* _data;

public:
	/**
	 * Creates a new <code>Heap</code> with the specified size.
	 * @param n	the maximum size of the heap.
	 */
	Heap(int n): _size(n), _nitem(0)
	{
		if (n > 0) { _data = new(Theap<T>[n]); }
		else { }
	}

	~Heap() { delete[] _data; }

	/**
	 * Checks whether the heap is empty.
	 * @return true if the heap empty.
	 */
	bool Empty() { return _nitem == 0; }

	/**
	 * Checks whether the heap is full.
	 * @return true if the heap full.
	 */
	bool Full() { return _nitem == _size - 1; }

	/**
	 * Returns the current number of elements if the heap.
	 * @return the current size of the heap
	 */
	int Nrank() { return _nitem; }

	/**
	 * Resets the heap (-> Nrank() = 0).
	 */
	void Reset() { _nitem = 0; }

	/**
	 * Pushes a new element in the heap with the specified key.
	 * @param item	the element to be inserted.
	 * @param ikey	the key that identify the element in the heap.
	 * @return true if the insertion is ok, fase if the heap is full.
	 */
	bool Push(const T& item, float ikey)
	{
		if (Full()) { return false; }
		int child = _nitem++;
		int parent = (child - 1) / 2;

		while (child > 0 && ikey < _data[parent].key) {
			_data[child] = _data[parent];
			child = parent;
			parent = (parent - 1) / 2;
		}
		_data[child].key = ikey;
		_data[child].element = item;

		return true;
	}

	/**
	 * Removes and returns the next element with the minimum key value.
	 * @param kitem	use to return the value of the key of the next element.
	 * @return	the element.
	 */
	T Pop(float* kitem = NULL)
	{
		if (Empty()) { return _data[0].element; }


		// Get top of the heap.
		T item = _data[0].element;
		if (kitem != NULL) { *kitem = _data[0].key; }

		--_nitem;
		int parent = 0;
		int child = 1;

		// Move up the lower branch of the heap.
		while (child < _nitem) {
			if (child + 1 < _nitem && _data[child].key > _data[child + 1].key) { child++; }

			if (_data[_nitem].key <= _data[child].key) { break; }

			// Move up child as parent.
			_data[parent] = _data[child];
			parent = child;
			child = child * 2 + 1;
		}
		_data[parent] = _data[_nitem];
		return item;
	}

	/**
	 * Removes from the heap the first occurence of the element
	 * specified by its key.
	 * @param kitem	the key of the element to be removed.
	 * @return	the element.
	 */
	Errc Remove(float kitem)
	{
		bool result = false;
		int i = 0;

		for (; i < _nitem && _data[i].key != kitem; i++) { ; }
		if (i >= _nitem) { return false; }

		--_nitem;
		int parent = i;
		int child = i * 2 + 1;

		// Move up the lower branch of the heap.
		while (child < _nitem) {
			if (child + 1 < _nitem && _data[child].key > _data[child + 1].key) { child++; }

			if (_data[_nitem].key <= _data[child].key) { break; }

			// Move up child as parent.
			_data[parent] = _data[child];
			parent = child;
			child = child * 2 + 1;
		}
		_data[parent] = _data[_nitem];

		return true;
	}
};

#endif // __PPLISTH__
