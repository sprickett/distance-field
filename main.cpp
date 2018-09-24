#include "distance_field.hpp"

int main(int argc, char* argv[])
{
    enum{H = 16, F = H*2, R2 = 196};
    DistanceFieldGenerator gen;

    std::vector<SrcType> data(F*F, 0);
    for(int y=0;y<F;++y)
        for(int x=0;x<F;++x)
        {
            int dx = x - H;
            int dy = y - H;
            data[x + y*F] = dx*dx + dy*dy < R2;
        }
            

    RegionOfInterest<SrcType> roi;
    roi.data = data.data();
    roi.width = F;
    roi.height = F;
    roi.stride = F;

    gen(roi,roi,0);

    return 0;
}