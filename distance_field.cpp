#include "distance_field.hpp"

PositionOrigin k_lut[] = {
    // # 21
    //   84
    
    // 0 00 --
    //   00 --

    // 1 01 11
    //   00 21
    { { 0, 0}, { 1, 0} },
    { { 0, 0}, { 0, 1} },
    { { 0, 1}, { 0,-1} },
    { { 0, 0}, { 0, 1} },

    // 2 10 11
    //   00 12
    
    // 3 11 11
    //   00 11

    // 4 00 21
    //   01 11

    // 5 01 11
    //   01 11

    // 6 10 11
    //   01 11

    // 7 11 12
    //   01 11

    // 8 00 12
    //   10 11

    // 9 01 11
    //   10 11

    // A 10 11
    //   10 11

    // B 11 21
    //   10 11

    // C 00 11
    //   11 11

    // D 01 11
    //   11 12

    // E 10 11
    //   11 21

    // F 11 --
    //   11 --
};

    Point kTbl[16][4] = {
        {   {0,0}, {0,0}, 
            {0,0}, {0,0}    },
        {   {0,0}, {0,0}, 
            {0,0}, {0,0}    },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
        { 
            {0,0}, {0,0}, 
            {0,0}, {0,0} 
        },
    };

    template<typename T>
    void setZero(T* data, DimType step, DistType count)
    {
        for(const T* end = data + step * count; data!=end;data+=step)
            *data = 0;
    }

    int DistanceFieldGenerator::operator()(const RegionOfInterest<SrcType>& src, RegionOfInterest<SrcType>& tgt, SrcType threshold )
    {
        if(src.width != tgt.width|| src.height != tgt.height)
            return -1;

        // reallocate cost memory
        unsigned w = src.width+1;
        unsigned h = src.height+2;
        cost_.clear();
        cost_.resize(w*h, max_l2_);
        RegionOfInterest<DistType> cost;
        cost.data = &cost_[w];
        cost.width = src.width;
        cost.height = src.height;
        cost.stride = w;
        
        // set cost border to zero
        setZero( &cost_[0],1,w );
        setZero( &cost_[w*h-w], 1, w);
        setZero( &cost_[w-1], w, h-2);

        // find seeds 
        for(unsigned y=0;y<src.height-1;++y)
        {
            DistType* C = cost.data + y*cost.stride; 
            const SrcType* S = src.data + y*src.stride;
            unsigned K = (threshold < S[0]) | (threshold < S[src.stride]) * 4;

            for(unsigned x=1; x<src.width; ++x )
            {
                K = (K & 5) * 2 | (threshold < S[x]) | (threshold < S[src.stride + x]) * 4;
                if(0<K & K<15)
                {
                    printf("%x",K);
                }
                else
                {
                    printf(" ");
                }
            }
            printf("\n"); 
        }
        unsigned K = (threshold < src.data[0]) |  (threshold < src.data[src.stride]) * 4;



    }

