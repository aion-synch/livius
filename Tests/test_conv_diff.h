#ifndef TEST_CONV_DIFF_H
#define TEST_CONV_DIFF_H

namespace corenc
{
    class test_conv_diff
    {
    public:
        test_conv_diff(){};
        ~test_conv_diff(){};
        void conv_diff_fem(const int h_ref_max, const int p_ref_max = 1) const;
        void conv_diff_eigen(const int h_ref_max, const int p_ref_max = 1) const;
    };
}

#endif // TEST_CONV_DIFF_H
