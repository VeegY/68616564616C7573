#include "../src/include/slicedvector.hpp"

int slicedvector_test()
{
    size_t dim_global(1000);
    SlivedVector slvec<double>(dim_global);

    if (slvec.get_dim_global() != dim_global)
        return false;

    for (size_t i(0); i < dim_global; ++i)
        slvec.set_global(i, static_cast<double>(i));

    SlicedVector slvec2<double>(5);
    slvec2 = slvec;

    return 0;
}

int main()
{
    return slicedvector_test();
}
