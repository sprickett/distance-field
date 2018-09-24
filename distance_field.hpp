#pragma once
#include <vector>
#include <array>

typedef unsigned char SrcType;
typedef unsigned int DistType;
typedef short DimType;

struct Point
{
    DimType x,y;
};
struct PositionOrigin
{
    Point position, origin;
};
template<typename T> struct RegionOfInterest
{
    T* data;
    DimType width;
    DimType height;
    DimType stride;
};

class DistanceFieldGenerator
{
private:
public:

    DistanceFieldGenerator(void)
        :max_l2_(DistType(-1))
        ,width_(0)
        ,height_(0)
        ,stride_(0)
    {}

    int operator()(const RegionOfInterest<unsigned char>& src, RegionOfInterest<unsigned char>& tgt, unsigned char threshold );

private:

    DistType max_l2_;
    DimType width_;
    DimType height_;
    DimType stride_;

    std::vector<DistType> cost_;
    std::array<std::vector<PositionOrigin>, 2> wave_;


};