#include <vector>
#include <set>
#include <cmath>

//Point definition
typedef std::pair<double, double> Point;

inline double dist(const Point& p1, const Point& p2){
    return std::sqrt(
	    (p1.first - p2.first) * (p1.first - p2.first)
	    + (p1.second - p2.second)*(p1.second-p2.second)
	    );
}

//Used to track lines of zeros of Phi_x
struct Tracker{
    Point head;
    double grad_value;
    double max_reach;
    std::vector<Point> zeros;

    Tracker(const Point& p, const double val, const double closeness){
	max_reach = closeness;
	head = p;
	grad_value = val;
	zeros = std::vector<Point>{};
    }

    inline bool Add(const Point& p, const double val){
	//Tries to add a point to tracker
	//Returns true if successful, false otherwise 
	//(If point is too far)
	//Also checks if Phi_y has changed sign, if yes
	//Adds a zero to zeros list
	if(dist(p, head) > max_reach) return false;
	
	if(val*grad_value < 0.){
	    //New zero
	    auto new_zero = Point(0.5*(p.first +head.first), 0.5*(p.second+head.second));
	    zeros.emplace_back(new_zero);
	}
	head = p;
	grad_value = val;
	return true;
    }
};

inline bool operator<(const Tracker& t1, const Tracker& t2){
    return t1.head < t2.head;
}
typedef std::vector<Tracker> TrackerList;

inline void TrackPoint(TrackerList& trackers, const Point& p, const double val, const double closeness){
    auto added = false;
    //Try to pt to exisiting tracker
    for(auto& tracker : trackers){
	if(tracker.Add(p, val)){
	    added = true;
	    break;
	}
    }
    if(added) return;

    //Otherwise create new tracker
    auto new_tracker = Tracker(p, val, closeness);
    trackers.emplace_back(new_tracker);
}

inline std::vector<Point> ExtractZeros(const TrackerList& trackers){
    auto all_zeros = std::vector<Point>{};
    for(const auto& tracker : trackers){
	for(const auto& pt : tracker.zeros)
	    all_zeros.emplace_back(pt);
    }
    return all_zeros;
}


