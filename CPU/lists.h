
#include "errc.h"
#include "generalParam.h"

#ifndef __PPLISTH__
#define __PPLISTH__

/*
 * Base class of complex data structures
 */
struct Slist {
public:  
   Slist *next;

   Slist(){next=NULL;}
};

/*
 * Base class of list.
 */
struct Blist {
protected:
   Slist *head;
   Slist *queue;

public:	
   Blist(){head=queue=NULL;}
   ~Blist(){}
   bool	Insert_queue(Slist* a){
      if (head!=NULL){
	 queue->next=a;
	 queue=a;
      }else{
	 head=queue=a;
      }
      return true;
   }

   bool	Insert_head(Slist* a){
      if (head!=NULL){
	 a->next=head;
	 head=a;
      }else
	 head=queue=a;
      return true;
   }
   Slist *Head(){return head;}
   Slist *Next(){return head;}
   Slist *Get_head(){
      if (head!=NULL){
	 Slist	*a;
	 a = head;
	 head = head->next;
	 if (head==NULL) queue=NULL;
	 return(a);
      }
      return(NULL);
   }
   bool Empty(){return(head==NULL);}
   bool Full(){return false;}
};

/*
 * Base class of items.
 */
template <class T>
struct Tlist : public Slist{
protected:
   T element;
public:
   Tlist(const T &a) : element(a) {}
   Tlist<T>* Next() {return (Tlist<T>*)next;}
   bool HasNext() {return next != NULL;}
   T Item() {return element;}
};

/**
 * The class <code>Fifo</code> defines a list of 
 * any primitive types (Uchar, Char, Slong, Float, ..)
 * that is a first in first out list.
 * Element is inserted by <code>push()</code> and
 * removed by <code>pop()</code>.
 */
template <class T>
class Fifo : public Blist {
public:

   /**
    * Gets without removing the next to be returned.
    * @return the next element to be poped.
    */
   Tlist<T> *Head(){ return (Tlist<T>*)Blist::Head(); }

   /**
    * Pushs the specified element in the Fifo.
    * Is is not a distinct copy of the element but a reference on it.
    * @param element	the element to be pushed.
    */
   void Push(const T &element){
      Blist::Insert_queue(new Tlist<T>(element));
   }

   /**
    * Gets and removes the nxt element in the Fifo (the first in).
    * @return the first element.
    */
   T Pop(){
      Tlist<T>* l=(Tlist<T>*)Blist::Get_head();
      T a= l->Item();
      delete l;
      return a;
   }

   /**
    * Resets the content of the fifo.
    */
   void Reset() { while (!Empty()) (void)Pop(); } 
};

/**
 * The class <code>Fifo</code> defines a list of 
 * any primitive types (Uchar, Char, Slong, Float, ..)
 * that is a last in first out list.
 * Element is inserted by <code>push()</code> and
 * removed by <code>pop()</code>.
 */
template <class T>
class Lifo :public Blist{
public :

   /**
    * Gets without removing the next to be returned.
    * @return the next element to be poped.
    */
   Tlist<T> *Head(){ return (Tlist<T>*)Blist::Head(); }

   /**
    * Pushs the specified element in the Lifo.
    * Is is not a distinct copy of the element but a reference on it.
    * @param element	the element to be pushed.
    */
   void Push(const T &element){
      Blist::Insert_head(new Tlist<T>(element));
   }

   /**
    * Gets and removes the next element in the Lifo (the first in).
    * @return the last element.
    */
   T Pop(){
      Tlist<T>* l=(Tlist<T>*)Blist::Get_head();
      T a= l->Item();
      delete l;
      return a;
   }

   /**
    * Resets the content of the lifo.
    */
   void Reset() { while (!Empty()) (void)Pop(); } 
};

/*
 * Base class of Heap.
 */
template <class T>
class Theap{
public :
    float key;
    T element;
	
    Theap(){}
    ~Theap(){}
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
class Heap{
private:
   int _size;
   int _nitem;
   Theap<T> *_data;

public:

   /**
    * Creates a new <code>Heap</code> with the specified size.
    * @param n	the maximum size of the heap.
    */
   Heap(int n):_size(n),_nitem(0) {
      if (n>0) _data=new(Theap<T>[n]);
      else{

	  }
   }

   ~Heap(){ delete[] _data; }

   /**
    * Checks whether the heap is empty.
    * @return true if the heap empty.
    */
   bool Empty(){ return(_nitem==0); }

   /**
    * Checks whether the heap is full.
    * @return true if the heap full.
    */
   bool Full(){return(_nitem==_size-1);}

   /**
    * Returns the current number of elements if the heap.
    * @return the current size of the heap
    */
   int Nrank(){return(_nitem);}

   /**
    * Resets the heap (-> Nrank() = 0).
    */
   void Reset() { _nitem=0;} 

   /**
    * Pushes a new element in the heap with the specified key.
    * @param item	the element to be inserted.
    * @param ikey	the key that identify the element in the heap.
    * @return true if the insertion is ok, fase if the heap is full.
    */
	bool Push(const T &item, float ikey){
		if (Full()){

			return false;
		}
		int child=_nitem++;
		int parent=(child-1)/2;

		while (child>0 && ikey<_data[parent].key){
			_data[child]=_data[parent];
			child=parent;
			parent=(parent-1)/2;
		}
		_data[child].key=ikey;
		_data[child].element=item;

		return true;
	} 

   /**
    * Removes and returns the next element with the minimum key value.
    * @param kitem	use to return the value of the key of the next element.
    * @return	the element.
    */
	T Pop(float *kitem=NULL){
		if (Empty()){

			return(_data[0].element);
		}


		// Get top of the heap.
		T item=_data[0].element;
		if (kitem!=NULL) *kitem=_data[0].key;

		--_nitem;
		int parent=0;
		int child=1;

		// Move up the lower branch of the heap.
		while (child<_nitem){
			if (child+1<_nitem && _data[child].key>_data[child+1].key)
			child++;

			if (_data[_nitem].key<=_data[child].key)
			break;

			// Move up child as parent.
			_data[parent]=_data[child];
			parent=child;
			child=child*2+1;
		}
		_data[parent]=_data[_nitem];
		return item;
	}

   /**
    * Removes from the heap the first occurence of the element
    * specified by its key.
    * @param kitem	the key of the element to be removed.
    * @return	the element.
    */
   Errc Remove(float kitem) {
      bool result=false;
      int i=0;

      for (; i<_nitem && _data[i].key!=kitem; i++);
      if (i>= _nitem) return false;

      --_nitem;
      int parent=i;
      int child=i*2+1;

      // Move up the lower branch of the heap.
      while (child<_nitem){
	 if (child+1<_nitem && _data[child].key>_data[child+1].key)
	    child++;

	 if (_data[_nitem].key<=_data[child].key)
	    break;

	 // Move up child as parent.
	 _data[parent]=_data[child];
	 parent=child;
	 child=child*2+1;
      }
      _data[parent]=_data[_nitem];

      return true;
   }
};

/**
 * The class <code>OrderedFifo</code> defines a ordered
 * list of element, which is implemented as an array of fifo.
 * Each element is indexed by a key value (a Long value).
 * Elements with the same key value are ordered in a fifo.
 * To create the new Heap named "foo" with a maximum of 100 floats
 * use: <code>Heap<Float> foo(100)</code>.
 */
template< typename T >
class OrderedFifo {
private:
   Long nfifo;
   Long nitem;
   Long min;
   Long max;

   Fifo< T > *tfifo;

public:
   
   /**
    * Creates a new ordered list with the specified maximum size
    * that corresponds to the maximum key value for elements.
    * @param size	the number of fifo.
    */
   OrderedFifo(Long size)/*: nfifo(t), nitem(0), min(0),max(0)*/ {
	   nfifo=t;
	   nitem=0;
	   min=0;
	   max=0;
      if (t) {
	 tfifo=(Fifo< T >*)calloc(size,sizeof(Fifo< T >));
	 if (tfifo==NULL) {

	 }
      }else
	 tfifo=NULL;
   }

   ~OrderedFifo() {
      for (Long i=0;i<nfifo;i++)
	 while (!(tfifo[i].Empty()))
	    tfifo[i].Pop();
      free(tfifo);
   }

   /**
    * Returns the current number of elements if the heap.
    * @return the current size of the heap
    */
   int Nrank(){ return nitem; }

   /**
    * Returns the i th fifo in the list. 
    * @param i	the index of the fifo.
    * @return	the related fifo or NULL.
    */
   Fifo< T > operator[](Long i) { return(tfifo[i]); }

   /**
    * Returns the index of the first not empty fifo.
    * @return the index in the list.  
    */
   Long First() {
      Long i;
      for (i=0;(i<nfifo) && (tfifo[i].Empty());i++);
      return(i);
   }
   
   /**
    * Returns the index of the last not empty fifo.
    * @return the index in the list.  
    */
   Long Last() {
      Long i;
      
      if (!nfifo) return(0);
      for (i=nfifo-1;(i>=0) && (tfifo[i].Empty());i--);
      return(i);
   }

   /**
    * Pushes a new element in the heap with the specified key.
    * @param item	the element to be inserted.
    * @param kitem	the key that identify the element in the ordered list.
    * @return true if the insertion is ok, fase if the heap is full.
    */
   bool Push(const T &item,Long kitem) {
      if (kitem>=nfifo)
	 return false;
      tfifo[kitem].Push(item);
      if (((min==max) && (min==0)) || (min>kitem)) min=kitem;
      if (max<kitem) max=kitem;
      nitem++;
      return true;
   }
   
   /**
    * Removes and returns the next element with the minimum key value.
    * @param kitem	use to return the value of the key of the next element.
    * @return	the related element.
    */
   T Pop(long *kitem=NULL) {
      T pt;
      if (Empty()) {

	 return T();
      }
      pt=tfifo[min].Pop();
      if (kitem!=NULL)
	 *kitem=min;
      if (tfifo[min].Empty()) {
	 if (max==min) min=max=0;
	 else min=First();
      }
      nitem--;
      return(pt);
   }
   
   /**
    * Checks whether the OrderedFifo is empty.
    * @return true if the OrderedFifo is empty, false otherwise.
    */
   bool Empty() { return((min==max) && (min==0)); }
};

#endif // __PPLISTH__
