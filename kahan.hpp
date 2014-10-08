#ifndef _KAHAN_HPP_
#define _KAHAN_HPP_

template<typename T>
struct kahan_sum
{
  T s,c,y,t;
  kahan_sum() : s(0.),c(0.),y(0.),t(0.){}
  
  kahan_sum<T> & operator=(const kahan_sum<T> & rhs)
  {
    this->s = rhs.s;
    this->c = rhs.c;
    this->y = rhs.y;
    this->t = rhs.t;
    return *this;
  }
  
  T & operator()( const T & v,const T & i )
  //! For use in algorithms (std::accumulate)
  {
    y=i-c;
    t=s+y;
    c=(t-s)-y;
    s=t;
    return s;
  }
  
  kahan_sum<T> & operator+=(const T & i)
  //!For accurate summation in blocks of code
  {
    y=i-c;
    t=s+y;
    c=(t-s)-y;
    s=t;
    return *this;
  }
  
  kahan_sum<T> & operator/=(const T & i) 
  {
    s /= i;
    return *this;
  }
  
  kahan_sum<T> operator/(const T & i) const
  {
    kahan_sum<T> temp;
    temp.s = this->s/i;
    return temp;
  }
  
  operator double() const
  {
    return s;
  }
};

#endif /* _KAHAN_HPP_ */
