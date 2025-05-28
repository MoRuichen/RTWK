#ifndef INTERVAL_H
#define INTERVAL_H

class interval {
  public:
  double min, max;
  
  interval() : min(+infinity),max(-infinity){}

  interval(double min,double max) :min(min),max(max){}

  double size() const{
    return max-min;
  }
  bool contain(double x) const{
    return x>=min && x<=max;
  }
  bool surrouds(double x) const{
    return x>min && x<max;
  }

  double clamp(double x) const{
    if(x<min) return min;
    if(x>max) return max;
    return x;
  }
  static const interval empty, universe;//static：表示这个成员属于类本身，而不是类的某个对象。比如所有快递柜（interval）共享同一个"空区间"的定义，而不是每个柜子自己存一份。意思就是这样写就说明有一个公用的概念，而不是每生成一个对象都生成一个无穷的值。
};


const interval interval::empty = interval (+infinity,-infinity);
const interval interval::universe = interval (-infinity,+infinity);
#endif