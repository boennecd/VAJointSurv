#pragma once

//  Concurrent queue of chapter 3, 
//  Used in the thread pool

#include <queue>
#include <mutex>

namespace cfaad {

template <class T>
class ConcurrentQueue
{
    template<class U>
    using lock_guard = std::lock_guard<U>;
    
    
    std::queue<T> myQueue;
	mutable std::mutex myMutex;
	std::condition_variable myCV;
	bool myInterrupt;

public:

	ConcurrentQueue() : myInterrupt(false) {}
	~ConcurrentQueue() { interrupt(); }

	bool empty() const
	{
		//	Lock
		lock_guard<std::mutex> lk(myMutex);
		//	Access underlying queue
		return myQueue.empty();
	}	//	Unlock

    //	Pop into argument
	bool tryPop(T& t)
	{
		//	Lock
		lock_guard<std::mutex> lk(myMutex);
		if (myQueue.empty()) return false;
		//	Move from queue
		t = move(myQueue.front());
		//	Combine front/pop
		myQueue.pop();

		return true;
	}	//	Unlock

    //	Pass t byVal or move with push( move( t))
	void push(T t)
	{
		{
			//	Lock
			lock_guard<std::mutex> lk(myMutex);
			//	Move into queue
			myQueue.push(move(t));
		}	//	Unlock before notification

        //	Unlock before notification 
		myCV.notify_one();
	}

	//	Wait if empty
	bool pop(T& t)
	{
		//	(Unique) lock
		std::unique_lock<std::mutex> lk(myMutex);

		//	Wait if empty, release lock until notified 
		while (!myInterrupt && myQueue.empty()) myCV.wait(lk);

		//	Re-acquire lock, resume 

		//  Check for interruption
		if (myInterrupt) return false;

		//	Combine front/pop 
		t = move(myQueue.front());
		myQueue.pop();

		return true;

	}	//	Unlock

	void interrupt()
	{
        {
            lock_guard<std::mutex> lk(myMutex);
            myInterrupt = true;
        }
		myCV.notify_all();
	}

    void resetInterrupt()
    {
        myInterrupt = false;
    }

    void clear()
    {
        std::queue<T> empty;
        swap(myQueue, empty);
    }
};

} // namespace cfaad