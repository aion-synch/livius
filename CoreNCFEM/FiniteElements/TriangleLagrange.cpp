#include "Triangle.h"
#include <iostream>
#include <algorithm>
#include <random>
#include "../../CoreNCA/Matrix.h"
#include "../../CoreNCA/MatrixSkyline.h"
using namespace corenc;
using namespace Mesh;
using namespace std;
using namespace Algebra;

CTriangleLagrangeBasis::CTriangleLagrangeBasis()
{
    m_det = 0;
    m_alpha[0][0] = m_alpha[0][1] = m_alpha[0][2] =
    m_alpha[1][0] = m_alpha[1][1] = m_alpha[1][2] =
    m_alpha[2][0] = m_alpha[2][1] = m_alpha[2][2] = 0;
    m_sp = 0;
    m_s = 0;
    m_1dorder = 0;
    m_order = 0;
    createS();
    m_normal = Point(0, 0, 0);
    m_number = 0;
    //m_all.resize(1);
}
CTriangleLagrangeBasis::CTriangleLagrangeBasis(const Point& p1, const Point& p2, const Point& p3, const int order)
{
    compD(p1, p2, p3);
    compAlpha(p1, p2, p3);
    compNormal(p1, p2, p3);
    m_order = order;
    createS();
    m_1dorder = 1;
    m_number = 3;
    //m_all.resize(1);
}

CTriangleLagrangeBasis::CTriangleLagrangeBasis(const Point* p, const int order)
{
    compD(p[0], p[1], p[2]);
    compAlpha(p[0], p[1], p[2]);
    compNormal(p[0], p[1], p[2]);
    m_order = order;
    createS();
    m_1dorder = 1;


}

const int CTriangleLagrangeBasis::createS()
{
    if (m_order < 2)
    {
        m_s = 3;
        m_sp = 0;
        return 0;
    }
    m_sp = std::max(int(3 * (m_order - 1) + (m_order - 2) * (m_order - 3) / 2), 0);
    m_s = 3 * (m_order) + (m_order - 1) * (m_order - 2) / 2;
    return 0;
}

const int CTriangleLagrangeBasis::GetNumberOfShapeFunctions() const
{
    return m_s;
}

const Point CTriangleLagrangeBasis::GetNormal() const
{
    return m_normal;
}

void CTriangleLagrangeBasis::ReverseNormal()
{
    m_normal.x = -m_normal.x;
    m_normal.y = -m_normal.y;
    m_normal.z = -m_normal.z;
}
const double CTriangleLagrangeBasis::m_L(const int k, const Point &p) const
{
    return (m_alpha[k][0] + m_alpha[k][1] * p.x + m_alpha[k][2] * p.y);
}

const double CTriangleLagrangeBasis::m_xi(const int k, const Point &p) const
{
    return std::pow(p.x, (k - 2));
}

const double CTriangleLagrangeBasis::GetShapeFunction(const int k, const Point &p) const
{
    switch (m_order)
    {
        case 1:
            return m_L(k, p);
        case 2:
            if (k < 3)
                return m_L(k, p) * (2 * m_L(k, p) - 1);
            switch (k)
            {
                case 3: return 4 * m_L(0, p) * m_L(1, p);
                case 4: return 4 * m_L(1, p) * m_L(2, p);
                case 5: return 4 * m_L(0, p) * m_L(2, p);
            }
        case 3:
            if (k < 3)
                return 0.5 * m_L(k, p) * (3 * m_L(k, p) - 1) * (3 * m_L(k, p) - 2);
            switch (k)
            {
                case 3: return 9. / 2. * m_L(0, p) * m_L(1, p) * (3 * m_L(0, p) - 1);
                case 6: return 9. / 2. * m_L(0, p) * m_L(1, p) * (3 * m_L(1, p) - 1);
                case 4: return 9. / 2. * m_L(1, p) * m_L(2, p) * (3 * m_L(1, p) - 1);
                case 7: return 9. / 2. * m_L(1, p) * m_L(2, p) * (3 * m_L(2, p) - 1);
                case 8: return 9. / 2. * m_L(2, p) * m_L(0, p) * (3 * m_L(2, p) - 1);
                case 5: return 9. / 2. * m_L(2, p) * m_L(0, p) * (3 * m_L(0, p) - 1);
                case 9: return 27. * m_L(0, p) * m_L(1, p) * m_L(2, p);
            }

        case 4:
            if (k < 3)
                return (4 * m_L(k, p) - 3) * (4 * m_L(k, p) - 1) * (2 * m_L(k, p) - 1) * m_L(k, p) / 3;
            switch (k)
            {
                case 6: return 4 * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 1) * (4 * m_L(1, p) - 1);
                case 7: return 4 * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 1) * (4 * m_L(2, p) - 1);
                case 8: return 4 * m_L(2, p) * m_L(0, p) * (4 * m_L(2, p) - 1) * (4 * m_L(0, p) - 1);
                /*case 6: return (8 * m_L(0, p) * m_L(0, p) - 2 * m_L(0, p)) * (8 * m_L(1, p) * m_L(1, p) - 2 * m_L(1, p));
                case 7: return (8 * m_L(1, p) * m_L(1, p) - 2 * m_L(1, p)) * (8 * m_L(2, p) * m_L(2, p) - 2 * m_L(2, p));
                case 8: return (8 * m_L(0, p) * m_L(0, p) - 2 * m_L(0, p)) * (8 * m_L(2, p) * m_L(2, p) - 2 * m_L(2, p));*/

                case 3: return 8. / 3. * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) * m_L(1, p);
                case 9: return 8. / 3. * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) * m_L(0, p);

                case 4: return 8. / 3. * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) * m_L(2, p);
                case 10: return 8. / 3. * m_L(2, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) * m_L(1, p);

                case 11: return 8. / 3. * m_L(2, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) * m_L(0, p);
                case 5: return 8. / 3. * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) * m_L(2, p);

                case 12: return 32 * m_L(0, p) * (4 * m_L(0, p) - 1) * m_L(1, p) * m_L(2, p);
                case 13: return 32 * m_L(1, p) * (4 * m_L(1, p) - 1) * m_L(0, p) * m_L(2, p);
                case 14: return 32 * m_L(2, p) * (4 * m_L(2, p) - 1) * m_L(1, p) * m_L(0, p);
                default:
                    break;
            }
        }
    return 0.;
}

const Point CTriangleLagrangeBasis::GetGradShapeFunction(const int k, const Point & p) const
{
    double temp = 0.;
    switch (m_order)
    {
        case 1:
            return Point(m_alpha[k][1], m_alpha[k][2]);
        case 2:
            if (k < 3)
                return Point{4 * m_alpha[k][1] * m_L(k, p) - m_alpha[k][1], 4 * m_alpha[k][2] * m_L(k, p) - m_alpha[k][2]};
            switch (k)
            {
                case 3: return Point{4 * (m_alpha[0][1] * m_L(1, p) + m_alpha[1][1] * m_L(0, p)), 4 * (m_alpha[0][2] * m_L(1, p) + m_alpha[1][2] * m_L(0, p))};
                case 4: return Point{4 * (m_alpha[2][1] * m_L(1, p) + m_alpha[1][1] * m_L(2, p)), 4 * (m_alpha[2][2] * m_L(1, p) + m_alpha[1][2] * m_L(2, p))};
                case 5: return Point{4 * (m_alpha[0][1] * m_L(2, p) + m_alpha[2][1] * m_L(0, p)), 4 * (m_alpha[0][2] * m_L(2, p) + m_alpha[2][2] * m_L(0, p))};
            }
        case 3:
            if (k < 3)
            {
                temp = m_L(k, p) * (27. * m_L(k, p) - 18.) / 2.;
                return Point(m_alpha[k][1] * (temp + 1), m_alpha[k][2] * (temp + 1));
            }
            switch (k)
            {
                case 3:
                    return Point(9. / 2. * (m_alpha[0][1] * m_L(1, p) * (3 * m_L(0, p) - 1) + m_alpha[1][1] * m_L(0, p) * (3 * m_L(0, p) - 1) + 3 * m_alpha[0][1] * m_L(0, p) * m_L(1, p)),
                            9. / 2. * (m_alpha[0][2] * m_L(1, p) * (3 * m_L(0, p) - 1) + m_alpha[1][2] * m_L(0, p) * (3 * m_L(0, p) - 1) + 3 * m_alpha[0][2] * m_L(0, p) * m_L(1, p)));
                case 6:
                    return Point(9. / 2. * (m_alpha[0][1] * m_L(1, p) * (3 * m_L(1, p) - 1) + m_alpha[1][1] * m_L(0, p) * (3 * m_L(1, p) - 1) + 3 * m_alpha[1][1] * m_L(0, p) * m_L(1, p)),
                            9. / 2. * (m_alpha[0][2] * m_L(1, p) * (3 * m_L(1, p) - 1) + m_alpha[1][2] * m_L(0, p) * (3 * m_L(1, p) - 1) + 3 * m_alpha[1][2] * m_L(0, p) * m_L(1, p)));
                case 4:
                    return Point(9. / 2. * (m_alpha[1][1] * m_L(2, p) * (3 * m_L(1, p) - 1) + m_alpha[2][1] * m_L(1, p) * (3 * m_L(1, p) - 1) + 3 * m_alpha[1][1] * m_L(1, p) * m_L(2, p)),
                            9. / 2. * (m_alpha[1][2] * m_L(2, p) * (3 * m_L(1, p) - 1) + m_alpha[2][2] * m_L(1, p) * (3 * m_L(1, p) - 1) + 3 * m_alpha[1][2] * m_L(1, p) * m_L(2, p)));
                case 7:
                    return Point(9. / 2. * (m_alpha[1][1] * m_L(2, p) * (3 * m_L(2, p) - 1) + m_alpha[2][1] * m_L(1, p) * (3 * m_L(2, p) - 1) + 3 * m_alpha[2][1] * m_L(1, p) * m_L(2, p)),
                            9. / 2. * (m_alpha[1][2] * m_L(2, p) * (3 * m_L(2, p) - 1) + m_alpha[2][2] * m_L(1, p) * (3 * m_L(2, p) - 1) + 3 * m_alpha[2][2] * m_L(1, p) * m_L(2, p)));
                case 8:
                    return Point(9. / 2. * (m_alpha[2][1] * m_L(0, p) * (3 * m_L(2, p) - 1) + m_alpha[0][1] * m_L(2, p) * (3 * m_L(2, p) - 1) + 3 * m_alpha[2][1] * m_L(0, p) * m_L(2, p)),
                            9. / 2. * (m_alpha[2][2] * m_L(0, p) * (3 * m_L(2, p) - 1) + m_alpha[0][2] * m_L(2, p) * (3 * m_L(2, p) - 1) + 3 * m_alpha[2][2] * m_L(0, p) * m_L(2, p)));
                case 5:
                    return Point(9. / 2. * (m_alpha[2][1] * m_L(0, p) * (3 * m_L(0, p) - 1) + m_alpha[0][1] * m_L(2, p) * (3 * m_L(0, p) - 1) + 3 * m_alpha[0][1] * m_L(0, p) * m_L(2, p)),
                            9. / 2. * (m_alpha[2][2] * m_L(0, p) * (3 * m_L(0, p) - 1) + m_alpha[0][2] * m_L(2, p) * (3 * m_L(0, p) - 1) + 3 * m_alpha[0][2] * m_L(0, p) * m_L(2, p)));
                case 9: return Point(27 * (m_alpha[0][1] * m_L(1, p) * m_L(2, p) + m_alpha[1][1] * m_L(0, p) * m_L(2, p) + m_alpha[2][1] * m_L(0, p) * m_L(1, p)), 27 * (m_alpha[0][2] * m_L(1, p) * m_L(2, p) + m_alpha[1][2] * m_L(0, p) * m_L(2, p) + m_alpha[2][2] * m_L(0, p) * m_L(1, p)));
                /*case 3: return Point(9. / 2. * (6 * m_alpha[0][1] * m_L(0, p) * m_L(1, p) + 3 * m_alpha[1][1] * m_L(0, p) * m_L(0, p) - m_alpha[0][1] * m_L(1, p) - m_alpha[1][1] * m_L(0, p)),
                        9. / 2. * (6 * m_alpha[0][2] * m_L(0, p) * m_L(1, p) + 3 * m_alpha[1][2] * m_L(0, p) * m_L(0, p) - m_alpha[0][2] * m_L(1, p) - m_alpha[1][2] * m_L(0, p)));
                case 4: return Point(9. / 2. * (6 * m_alpha[1][1] * m_L(1, p) * m_L(0, p) + 3 * m_alpha[0][1] * m_L(1, p) * m_L(1, p) - m_alpha[1][1] * m_L(0, p) - m_alpha[0][1] * m_L(1, p)),
                        9. / 2. * (6 * m_alpha[1][2] * m_L(1, p) * m_L(0, p) + 3 * m_alpha[0][2] * m_L(1, p) * m_L(1, p) - m_alpha[1][2] * m_L(0, p) - m_alpha[0][2] * m_L(1, p)));

                case 5: return Point(9. / 2. * (6 * m_alpha[1][1] * m_L(1, p) * m_L(2, p) + 3 * m_alpha[2][1] * m_L(1, p) * m_L(1, p) - m_alpha[1][1] * m_L(2, p) - m_alpha[2][1] * m_L(1, p)),
                        9. / 2. * (6 * m_alpha[1][2] * m_L(1, p) * m_L(2, p) + 3 * m_alpha[2][2] * m_L(1, p) * m_L(1, p) - m_alpha[1][2] * m_L(2, p) - m_alpha[2][2] * m_L(1, p)));
                case 6: return Point(9. / 2. * (6 * m_alpha[2][1] * m_L(2, p) * m_L(1, p) + 3 * m_alpha[1][1] * m_L(2, p) * m_L(2, p) - m_alpha[2][1] * m_L(1, p) - m_alpha[1][1] * m_L(2, p)),
                        9. / 2. * (6 * m_alpha[2][2] * m_L(2, p) * m_L(1, p) + 3 * m_alpha[1][2] * m_L(2, p) * m_L(2, p) - m_alpha[2][2] * m_L(1, p) - m_alpha[1][2] * m_L(2, p)));

                case 7: return Point(9. / 2. * (6 * m_alpha[2][1] * m_L(2, p) * m_L(0, p) + 3 * m_alpha[0][1] * m_L(2, p) * m_L(2, p) - m_alpha[2][1] * m_L(0, p) - m_alpha[0][1] * m_L(2, p)),
                        9. / 2. * (6 * m_alpha[2][2] * m_L(2, p) * m_L(0, p) + 3 * m_alpha[0][2] * m_L(2, p) * m_L(2, p) - m_alpha[2][2] * m_L(0, p) - m_alpha[0][2] * m_L(2, p)));
                case 8: return Point(9. / 2. * (6 * m_alpha[0][1] * m_L(0, p) * m_L(2, p) + 3 * m_alpha[2][1] * m_L(0, p) * m_L(0, p) - m_alpha[0][1] * m_L(2, p) - m_alpha[2][1] * m_L(0, p)),
                        9. / 2. * (6 * m_alpha[0][2] * m_L(0, p) * m_L(2, p) + 3 * m_alpha[2][2] * m_L(0, p) * m_L(0, p) - m_alpha[0][2] * m_L(2, p) - m_alpha[2][2] * m_L(0, p)));

                case 9: return Point(27 * (m_alpha[0][1] * m_L(1, p) * m_L(2, p) + m_alpha[1][1] * m_L(0, p) * m_L(2, p) + m_alpha[2][1] * m_L(0, p) * m_L(1, p)), 27 * (m_alpha[0][2] * m_L(1, p) * m_L(2, p) + m_alpha[1][2] * m_L(0, p) * m_L(2, p) + m_alpha[2][2] * m_L(0, p) * m_L(1, p)));*/
            }

        case 4:
            if(k < 3)
                //return Point(1. / 3 * (m_alpha[k][1] * (4. * m_L(k, p) - 3.) * (4. * m_L(k, p) - 1.) * (2. * m_L(k, p) - 1) + 4. * m_alpha[k][1] * m_L(k, p) * (4. * m_L(k, p) - 1.) * (2. * m_L(k, p) - 1) + 4 * m_alpha[k][1] * m_L(k, p) * (4. * m_L(k, p) - 3.) * (2. * m_L(k, p) - 1.) + 2. * m_alpha[k][1] * m_L(k, p) * (4. * m_L(k, p) - 3.) * (4. * m_L(k, p) - 1.)),
                  //          1. / 3 * (m_alpha[k][2] * (4. * m_L(k, p) - 3.) * (4. * m_L(k, p) - 1.) * (2. * m_L(k, p) - 1) + 4. * m_alpha[k][2] * m_L(k, p) * (4. * m_L(k, p) - 1.) * (2. * m_L(k, p) - 1) + 4 * m_alpha[k][2] * m_L(k, p) * (4. * m_L(k, p) - 3.) * (2. * m_L(k, p) - 1.) + 2. * m_alpha[k][2] * m_L(k, p) * (4. * m_L(k, p) - 3.) * (4. * m_L(k, p) - 1.)));
            {
                temp = (4 * m_L(k, p) - 3) * (4 * m_L(k, p) - 1) * (2 * m_L(k, p) - 1) + 4 * m_L(k, p) * (4 * m_L(k, p) - 1) * (2 * m_L(k, p) - 1) +
                        4 * m_L(k, p) * (4 * m_L(k, p) - 3) * (2 * m_L(k, p) - 1) + 2 * m_L(k, p) * (4 * m_L(k, p) - 3) * (4 * m_L(k, p) - 1);
                return Point(m_alpha[k][1] / 3 * temp, m_alpha[k][2] / 3 * temp);
            }
            if (k > 5 && k < 9)
            {
                int i, j;
                switch (k)
                {
                    case 6: i = 0; j = 1; break;
                    case 7: i = 1; j = 2; break;
                    case 8: i = 0; j = 2; break;
                }
                //return Point(4. * (m_alpha[i][1] * m_L(i, p) * (4. * m_L(i, p) - 1) * (4 * m_L(j, p) - 1) + m_alpha[j][1] * m_L(i, p) * (4 * m_L(i, p) - 1) * (4 * m_L(j, p) - 1) + 4 * m_alpha[i][1] * m_L(i, p) * m_L(j, p) * (4 - m_L(j, p) - 1) + 4 * m_alpha[j][1] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1)),
                //             4. * (m_alpha[i][2] * m_L(i, p) * (4. * m_L(i, p) - 1) * (4 * m_L(j, p) - 1) + m_alpha[j][2] * m_L(i, p) * (4 * m_L(i, p) - 1) * (4 * m_L(j, p) - 1) + 4 * m_alpha[i][2] * m_L(i, p) * m_L(j, p) * (4 - m_L(j, p) - 1) + 4 * m_alpha[j][2] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1)));
                return Point(4 * (m_alpha[i][1] * (4 * m_L(i, p) - 1) * m_L(j, p) * (4 * m_L(j, p) - 1) + m_alpha[j][1] * m_L(i, p) * (4 * m_L(i, p) - 1) * (4 * m_L(j, p) - 1) + 4 * m_alpha[i][1] * m_L(i, p) * m_L(j, p) * (4 * m_L(j, p) - 1) + 4 * m_alpha[j][1] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1)),
                             4 * (m_alpha[i][2] * (4 * m_L(i, p) - 1) * m_L(j, p) * (4 * m_L(j, p) - 1) + m_alpha[j][2] * m_L(i, p) * (4 * m_L(i, p) - 1) * (4 * m_L(j, p) - 1) + 4 * m_alpha[i][2] * m_L(i, p) * m_L(j, p) * (4 * m_L(j, p) - 1) + 4 * m_alpha[j][2] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1)));
            }
            if ((k > 2 && k < 6) || (k > 8 && k < 12))
            {
                int i, j;
                switch (k)
                {
                    case 3: i = 0; j = 1; break;
                    case 9: i = 1; j = 0; break;
                    case 4: i = 1; j = 2; break;
                    case 10: i = 2; j = 1; break;
                    case 11: i = 2; j = 0; break;
                    case 5: i = 0; j = 2; break;
                }
                return Point(8. / 3 * (m_alpha[i][1] * m_L(j, p) * (4 * m_L(i, p) - 1) * (4 * m_L(i, p) - 2) + m_alpha[j][1] * m_L(i, p) * (4 * m_L(i, p) - 1) * (4 * m_L(i, p) - 2) + 4 * m_alpha[i][1] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 2) + 4 * m_alpha[i][1] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1)),
                             8. / 3 * (m_alpha[i][2] * m_L(j, p) * (4 * m_L(i, p) - 1) * (4 * m_L(i, p) - 2) + m_alpha[j][2] * m_L(i, p) * (4 * m_L(i, p) - 1) * (4 * m_L(i, p) - 2) + 4 * m_alpha[i][2] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 2) + 4 * m_alpha[i][2] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1)));
            }
            int i, j, ij;
            switch (k)
            {
                case 12: i = 0; j = 1; ij = 2; break;
                case 13: i = 1; j = 0; ij = 2; break;
                case 14: i = 2; j = 1; ij = 0; break;
            }
            return Point(32. * (m_alpha[i][1] * m_L(j, p) * m_L(ij, p) * (4 * m_L(i, p) - 1) + m_alpha[j][1] * m_L(i, p) * m_L(ij, p) * (4 * m_L(i, p) - 1) + m_alpha[ij][1] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1) + 4 * m_alpha[i][1] * m_L(i, p) * m_L(j, p) * m_L(ij, p)),
                         32. * (m_alpha[i][2] * m_L(j, p) * m_L(ij, p) * (4 * m_L(i, p) - 1) + m_alpha[j][2] * m_L(i, p) * m_L(ij, p) * (4 * m_L(i, p) - 1) + m_alpha[ij][2] * m_L(i, p) * m_L(j, p) * (4 * m_L(i, p) - 1) + 4 * m_alpha[i][2] * m_L(i, p) * m_L(j, p) * m_L(ij, p)));
            switch (k)
            {
                /*case 6: return Point(2. * (m_alpha[0][1] * m_L(0, p) * (4. * m_L(0, p) - 1) * (4 * m_L(1, p) - 1) + m_alpha[1][1] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(1, p) - 1) + 4 * m_alpha[0][1] * m_L(0, p) * m_L(1, p) * (4 - m_L(1, p) - 1) + 4 * m_alpha[1][1] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 1)),
                        2. * (m_alpha[0][2] * m_L(0, p) * (4. * m_L(0, p) - 1) * (4 * m_L(1, p) - 1) + m_alpha[1][2] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(1, p) - 1) + 4 * m_alpha[0][2] * m_L(0, p) * m_L(1, p) * (4 - m_L(1, p) - 1) + 4 * m_alpha[1][2] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 1)));
                case 7: return Point(2. * (m_alpha[1][1] * m_L(1, p) * (4. * m_L(1, p) - 1) * (4 * m_L(2, p) - 1) + m_alpha[2][1] * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(2, p) - 1) + 4 * m_alpha[1][1] * m_L(1, p) * m_L(2, p) * (4 - m_L(2, p) - 1) + 4 * m_alpha[2][1] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 1)),
                    2. * (m_alpha[1][2] * m_L(1, p) * (4. * m_L(1, p) - 1) * (4 * m_L(2, p) - 1) + m_alpha[2][2] * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(2, p) - 1) + 4 * m_alpha[1][2] * m_L(1, p) * m_L(2, p) * (4 - m_L(2, p) - 1) + 4 * m_alpha[2][2] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 1)));
                case 8: return Point(2. * (m_alpha[0][1] * m_L(0, p) * (4. * m_L(0, p) - 1) * (4 * m_L(2, p) - 1) + m_alpha[2][1] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(2, p) - 1) + 4 * m_alpha[0][1] * m_L(0, p) * m_L(2, p) * (4 - m_L(2, p) - 1) + 4 * m_alpha[2][1] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 1)),
                    2. * (m_alpha[0][2] * m_L(0, p) * (4. * m_L(0, p) - 1) * (4 * m_L(2, p) - 1) + m_alpha[2][2] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(2, p) - 1) + 4 * m_alpha[0][2] * m_L(0, p) * m_L(2, p) * (4 - m_L(2, p) - 1) + 4 * m_alpha[2][2] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 1)));

                case 3: return Point(8. / 3 * (m_alpha[0][1] * m_L(1, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + m_alpha[1][1] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][1] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][1] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 1)),
                        8. / 3 * (m_alpha[0][2] * m_L(1, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + m_alpha[1][2] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][2] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][2] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 1)));
                case 9: return Point(8. / 3 * (m_alpha[1][1] * m_L(0, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + m_alpha[0][1] * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][1] * m_L(1, p) * m_L(0, p) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][1] * m_L(1, p) * m_L(0, p) * (4 * m_L(1, p) - 1)),
                        8. / 3 * (m_alpha[1][2] * m_L(0, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + m_alpha[0][2] * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][2] * m_L(1, p) * m_L(0, p) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][2] * m_L(1, p) * m_L(0, p) * (4 * m_L(1, p) - 1)));

                case 4: return Point(8. / 3 * (m_alpha[1][1] * m_L(2, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + m_alpha[2][1] * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][1] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][1] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 1)),
                        8. / 3 * (m_alpha[1][2] * m_L(2, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + m_alpha[2][2] * m_L(1, p) * (4 * m_L(1, p) - 1) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][2] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 2) + 4 * m_alpha[1][2] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 1)));
                case 10: return Point(8. / 3 * (m_alpha[2][1] * m_L(1, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + m_alpha[1][1] * m_L(2, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][1] * m_L(2, p) * m_L(1, p) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][1] * m_L(2, p) * m_L(1, p) * (4 * m_L(2, p) - 1)),
                    8. / 3 * (m_alpha[2][2] * m_L(1, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + m_alpha[1][2] * m_L(2, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][2] * m_L(2, p) * m_L(1, p) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][2] * m_L(2, p) * m_L(1, p) * (4 * m_L(2, p) - 1)));

                case 11: return Point(8. / 3 * (m_alpha[2][1] * m_L(0, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + m_alpha[0][1] * m_L(2, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][1] * m_L(2, p) * m_L(0, p) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][1] * m_L(2, p) * m_L(0, p) * (4 * m_L(2, p) - 1)),
                8. / 3 * (m_alpha[2][2] * m_L(0, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + m_alpha[0][2] * m_L(2, p) * (4 * m_L(2, p) - 1) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][2] * m_L(2, p) * m_L(0, p) * (4 * m_L(2, p) - 2) + 4 * m_alpha[2][2] * m_L(2, p) * m_L(0, p) * (4 * m_L(2, p) - 1)));
                case 5: return Point(8. / 3 * (m_alpha[0][1] * m_L(2, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + m_alpha[2][1] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][1] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][1] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 1)),
                    8. / 3 * (m_alpha[0][2] * m_L(2, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + m_alpha[2][2] * m_L(0, p) * (4 * m_L(0, p) - 1) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][2] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 2) + 4 * m_alpha[0][2] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 1)));*/

                /*case 12: return Point(32. * (m_alpha[0][1] * m_L(1, p) * m_L(2, p) * (4 * m_L(0, p) - 1) + m_alpha[1][1] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 1) + m_alpha[2][1] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 1) + 4 * m_alpha[0][1] * m_L(0, p) * m_L(1, p) * m_L(2, p)),
                        32. * (m_alpha[0][2] * m_L(1, p) * m_L(2, p) * (4 * m_L(0, p) - 1) + m_alpha[1][2] * m_L(0, p) * m_L(2, p) * (4 * m_L(0, p) - 1) + m_alpha[2][2] * m_L(0, p) * m_L(1, p) * (4 * m_L(0, p) - 1) + 4 * m_alpha[0][2] * m_L(0, p) * m_L(1, p) * m_L(2, p)));
                case 13: return Point(32. * (m_alpha[0][1] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 1) + m_alpha[1][1] * m_L(0, p) * m_L(2, p) * (4 * m_L(1, p) - 1) + m_alpha[2][1] * m_L(0, p) * m_L(1, p) * (4 * m_L(1, p) - 1) + 4 * m_alpha[1][1] * m_L(0, p) * m_L(1, p) * m_L(2, p)),
                    32. * (m_alpha[0][2] * m_L(1, p) * m_L(2, p) * (4 * m_L(1, p) - 1) + m_alpha[1][2] * m_L(0, p) * m_L(2, p) * (4 * m_L(1, p) - 1) + m_alpha[2][2] * m_L(0, p) * m_L(1, p) * (4 * m_L(1, p) - 1) + 4 * m_alpha[1][2] * m_L(0, p) * m_L(1, p) * m_L(2, p)));
                case 14: return Point(32. * (m_alpha[0][1] * m_L(1, p) * m_L(2, p) * (4 * m_L(2, p) - 1) + m_alpha[1][1] * m_L(0, p) * m_L(2, p) * (4 * m_L(2, p) - 1) + m_alpha[2][1] * m_L(0, p) * m_L(1, p) * (4 * m_L(2, p) - 1) + 4 * m_alpha[2][1] * m_L(0, p) * m_L(1, p) * m_L(2, p)),
                    32. * (m_alpha[0][2] * m_L(1, p) * m_L(2, p) * (4 * m_L(2, p) - 1) + m_alpha[1][2] * m_L(0, p) * m_L(2, p) * (4 * m_L(2, p) - 1) + m_alpha[2][2] * m_L(0, p) * m_L(1, p) * (4 * m_L(2, p) - 1) + 4 * m_alpha[2][2] * m_L(0, p) * m_L(1, p) * m_L(2, p)));*/

                /*case 0: temp = (256 * m_L(0, p) * m_L(0, p) * m_L(0, p) - 144 * m_L(0, p) * m_L(0, p) + 48 * m_L(0, p) - 3) / 3.;
                    return Point(temp * m_alpha[0][1], temp * m_alpha[0][2]);
                case 1: temp = (256 * m_L(1, p) * m_L(1, p) * m_L(1, p) - 144 * m_L(1, p) * m_L(1, p) + 48 * m_L(1, p) - 3) / 3.;
                    return Point(temp * m_alpha[0][1], temp * m_alpha[0][2]);
                case 2: temp = (256 * m_L(2, p) * m_L(2, p) * m_L(2, p) - 144 * m_L(2, p) * m_L(2, p) + 48 * m_L(2, p) - 3) / 3.;
                    return Point(temp * m_alpha[0][1], temp * m_alpha[0][2]);

                case 6: return Point(m_alpha[0][1] * (16 * m_L(0, p) - 2) * (8 * m_L(1, p) * m_L(1, p) - 2 * m_L(1, p)) +
                            m_alpha[1][1] * (16 * m_L(1, p) - 2) * (8 * m_L(0, p) * m_L(0, p) - 2 * m_L(0, p)),
                            m_alpha[0][2] * (16 * m_L(0, p) - 2) * (8 * m_L(1, p) * m_L(1, p) - 2 * m_L(1, p)) +
                                            m_alpha[1][2] * (16 * m_L(1, p) - 2) * (8 * m_L(0, p) * m_L(0, p) - 2 * m_L(0, p)));
                case 7: return Point(m_alpha[1][1] * (16 * m_L(1, p) - 2) * (8 * m_L(2, p) * m_L(2, p) - 2 * m_L(2, p)) +
                            m_alpha[2][1] * (16 * m_L(2, p) - 2) * (8 * m_L(1, p) * m_L(1, p) - 2 * m_L(1, p)),
                            m_alpha[1][2] * (16 * m_L(1, p) - 2) * (8 * m_L(2, p) * m_L(2, p) - 2 * m_L(2, p)) +
                                            m_alpha[2][2] * (16 * m_L(2, p) - 2) * (8 * m_L(1, p) * m_L(1, p) - 2 * m_L(1, p)));
                case 8: return Point(m_alpha[0][1] * (16 * m_L(0, p) - 2) * (8 * m_L(2, p) * m_L(2, p) - 2 * m_L(2, p)) +
                            m_alpha[2][1] * (16 * m_L(2, p) - 2) * (8 * m_L(0, p) * m_L(0, p) - 2 * m_L(0, p)),
                            m_alpha[0][2] * (16 * m_L(0, p) - 2) * (8 * m_L(2, p) * m_L(2, p) - 2 * m_L(2, p)) +
                                            m_alpha[2][2] * (16 * m_L(2, p) - 2) * (8 * m_L(0, p) * m_L(0, p) - 2 * m_L(0, p)));

                case 3: return Point(8. / 3 * (m_alpha[1][1] * (16 * m_L(0, p) * m_L(0, p) * m_L(0, p) - 12 * m_L(0, p) * m_L(0, p) + 2 * m_L(0, p)) + m_L(1, p) * m_alpha[0][1] * (48 * m_L(0, p) - 24 * m_L(0, p) + 2)),
                            8. / 3 * (m_alpha[1][2] * (16 * m_L(0, p) * m_L(0, p) * m_L(0, p) - 12 * m_L(0, p) * m_L(0, p) + 2 * m_L(0, p)) + m_L(1, p) * m_alpha[0][2] * (48 * m_L(0, p) - 24 * m_L(0, p) + 2)));
                case 9: return Point(8. / 3 * (m_alpha[0][1] * (16 * m_L(1, p) * m_L(1, p) * m_L(1, p) - 12 * m_L(1, p) * m_L(1, p) + 2 * m_L(1, p)) + m_L(0, p) * m_alpha[1][1] * (48 * m_L(1, p) - 24 * m_L(1, p) + 2)),
                            8. / 3 * (m_alpha[0][2] * (16 * m_L(1, p) * m_L(1, p) * m_L(1, p) - 12 * m_L(1, p) * m_L(1, p) + 2 * m_L(1, p)) + m_L(0, p) * m_alpha[1][2] * (48 * m_L(1, p) - 24 * m_L(1, p) + 2)));

                case 4: return Point(8. / 3 * (m_alpha[2][1] * (16 * m_L(1, p) * m_L(1, p) * m_L(1, p) - 12 * m_L(1, p) * m_L(1, p) + 2 * m_L(1, p)) + m_L(2, p) * m_alpha[1][1] * (48 * m_L(1, p) - 24 * m_L(1, p) + 2)),
                            8. / 3 * (m_alpha[2][2] * (16 * m_L(1, p) * m_L(1, p) * m_L(1, p) - 12 * m_L(1, p) * m_L(1, p) + 2 * m_L(1, p)) + m_L(2, p) * m_alpha[1][2] * (48 * m_L(1, p) - 24 * m_L(1, p) + 2)));
                case 10: return Point(8. / 3 * (m_alpha[1][1] * (16 * m_L(2, p) * m_L(2, p) * m_L(2, p) - 12 * m_L(2, p) * m_L(2, p) + 2 * m_L(2, p)) + m_L(1, p) * m_alpha[2][1] * (48 * m_L(2, p) - 24 * m_L(2, p) + 2)),
                            8. / 3 * (m_alpha[1][2] * (16 * m_L(2, p) * m_L(2, p) * m_L(2, p) - 12 * m_L(2, p) * m_L(2, p) + 2 * m_L(2, p)) + m_L(1, p) * m_alpha[2][2] * (48 * m_L(2, p) - 24 * m_L(2, p) + 2)));

                case 11: return Point(8. / 3 * (m_alpha[2][1] * (16 * m_L(0, p) * m_L(0, p) * m_L(0, p) - 12 * m_L(0, p) * m_L(0, p) + 2 * m_L(0, p)) + m_L(2, p) * m_alpha[0][1] * (48 * m_L(0, p) - 24 * m_L(0, p) + 2)),
                            8. / 3 * (m_alpha[2][2] * (16 * m_L(0, p) * m_L(0, p) * m_L(0, p) - 12 * m_L(0, p) * m_L(0, p) + 2 * m_L(0, p)) + m_L(2, p) * m_alpha[0][2] * (48 * m_L(0, p) - 24 * m_L(0, p) + 2)));
                case 5: return Point(8. / 3 * (m_alpha[2][1] * (16 * m_L(0, p) * m_L(0, p) * m_L(0, p) - 12 * m_L(0, p) * m_L(0, p) + 2 * m_L(0, p)) + m_L(2, p) * m_alpha[0][1] * (48 * m_L(0, p) - 24 * m_L(0, p) + 2)),
                            8. / 3 * (m_alpha[2][2] * (16 * m_L(0, p) * m_L(0, p) * m_L(0, p) - 12 * m_L(0, p) * m_L(0, p) + 2 * m_L(0, p)) + m_L(2, p) * m_alpha[0][2] * (48 * m_L(0, p) - 24 * m_L(0, p) + 2)));

                case 12: return Point(32 * (m_alpha[0][1] * m_L(1, p) * m_L(2, p) * (8 * m_L(0, p) - 1) + m_alpha[1][1] * m_L(2, p) * (4 * m_L(0, p) * m_L(0, p) - m_L(0, p)) + m_alpha[2][1] * m_L(1, p) * (4 * m_L(0, p) * m_L(0, p) - m_L(0, p))),
                            32 * (m_alpha[0][2] * m_L(1, p) * m_L(2, p) * (8 * m_L(0, p) - 1) + m_alpha[1][2] * m_L(2, p) * (4 * m_L(0, p) * m_L(0, p) - m_L(0, p)) + m_alpha[2][2] * m_L(1, p) * (4 * m_L(0, p) * m_L(0, p) - m_L(0, p))));
                case 13: return Point(32 * (m_alpha[1][1] * m_L(0, p) * m_L(2, p) * (8 * m_L(1, p) - 1) + m_alpha[0][1] * m_L(2, p) * (4 * m_L(1, p) * m_L(1, p) - m_L(1, p)) + m_alpha[2][1] * m_L(0, p) * (4 * m_L(1, p) * m_L(1, p) - m_L(1, p))),
                            32 * (m_alpha[1][2] * m_L(0, p) * m_L(2, p) * (8 * m_L(1, p) - 1) + m_alpha[0][2] * m_L(2, p) * (4 * m_L(1, p) * m_L(1, p) - m_L(1, p)) + m_alpha[2][2] * m_L(0, p) * (4 * m_L(1, p) * m_L(1, p) - m_L(1, p))));
                case 14: return Point(32 * (m_alpha[2][1] * m_L(0, p) * m_L(1, p) * (8 * m_L(2, p) - 1) + m_alpha[0][1] * m_L(1, p) * (4 * m_L(2, p) * m_L(2, p) - m_L(2, p)) + m_alpha[1][1] * m_L(0, p) * (4 * m_L(2, p) * m_L(2, p) - m_L(2, p))),
                            32 * (m_alpha[2][2] * m_L(0, p) * m_L(1, p) * (8 * m_L(2, p) - 1) + m_alpha[0][2] * m_L(1, p) * (4 * m_L(2, p) * m_L(2, p) - m_L(2, p)) + m_alpha[2][2] * m_L(0, p) * (4 * m_L(2, p) * m_L(2, p) - m_L(2, p))));*/
                default:
                    break;
        }
    }
    return Point(0, 0);
}

CTriangleLagrangeBasis::CTriangleLagrangeBasis(const CTriangleLagrangeBasis& t)
{
    m_normal = t.m_normal;
    m_det = t.m_det;
    m_alpha[0][0] = t.m_alpha[0][0];
    m_alpha[0][1] = t.m_alpha[0][1];
    m_alpha[0][2] = t.m_alpha[0][2];

    m_alpha[1][0] = t.m_alpha[1][0];
    m_alpha[1][1] = t.m_alpha[1][1];
    m_alpha[1][2] = t.m_alpha[1][2];

    m_alpha[2][0] = t.m_alpha[2][0];
    m_alpha[2][1] = t.m_alpha[2][1];
    m_alpha[2][2] = t.m_alpha[2][2];
    m_number = t.m_number;
    m_order = t.m_order;
    m_1dorder = t.m_1dorder;
    //m_w = t.m_w;
    m_s = t.m_s;
    m_sp = t.m_sp;
    m_all = t.m_all;
}

void CTriangleLagrangeBasis::compD(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
{
    double pz;
    //m_det = (p2.y - p1.y)*(p3.z - p1.z) - (p3.y - p1.y)*(p2.z - p1.z);
    //m_det *= m_det;
    //py = (p3.x - p1.x)*(p2.z - p1.z) - (p2.x - p1.x)*(p3.z - p1.z);
    //py *= py;
    pz = (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y);
    //pz *= pz;
    //m_det = sqrt(m_det + py + pz);
    m_det = pz;
}

void CTriangleLagrangeBasis::compAlpha(const Mesh::Point& p1, const Mesh::Point& p2, const Mesh::Point& p3)
{
    m_alpha[0][0] = (p2.x * p3.y - p2.y * p3.x) / m_det;
    m_alpha[1][0] = -(p1.x * p3.y - p1.y * p3.x) / m_det;
    m_alpha[2][0] = (p1.x * p2.y - p1.y * p2.x) / m_det;

    m_alpha[0][1] = -(p3.y - p2.y) / m_det;
    m_alpha[1][1] = (p3.y - p1.y) / m_det;
    m_alpha[2][1] = -(p2.y - p1.y) / m_det;

    m_alpha[0][2] = (p3.x - p2.x) / m_det;
    m_alpha[1][2] = -(p3.x - p1.x) / m_det;
    m_alpha[2][2] = (p2.x - p1.x) / m_det;
}

void CTriangleLagrangeBasis::compNormal(const Mesh::Point &p1, const Mesh::Point &p2, const Mesh::Point &p3)
{
    Point p;
    p.x = p1.y*(p2.z - p3.z) + p2.y*(p3.z - p1.z) + p3.y*(p1.z - p2.z);
    p.y = p1.z*(p2.x - p3.x) + p2.z*(p3.x - p1.x) + p3.z*(p1.x - p2.x);
    p.z = p1.x*(p2.y - p3.y) + p2.x*(p3.y - p1.y) + p3.x*(p1.y - p2.y);
    double t = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    p.x = (p.x) / t;
    p.y = (p.y) / t;
    p.z = (p.z) / t;
    m_normal = p;
}

const double CTriangleLagrangeBasis::GetValue(const Point& p) const
{
    return 0.;
}
const int CTriangleLagrangeBasis::IncreaseOrder()
{
    ++m_order;
    ++m_1dorder;
    createS();
    m_all.push_back(m_sp);
    return 0;
}
namespace wtf
{
    const Point mid_point(const Point& p1, const Point& p2)
    {
        return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
    }

    const Point s_point(const Point& p1, const Point& p2, const double s)
    {
        return Point((p1.x + p2.x) / s, (p1.y + p2.y) / s);
    }

    const Point center_point(const Point& p1, const Point& p2, const Point& p3)
    {
        return Point((p1.x + p2.x + p3.x) / 3, (p1.y + p2.y + p3.y) / 3);
    }
}
const double CTriangleLagrangeBasis::GetWeight(const int node, const vector<Point>& pts, const std::function<const double(const Point&)>& func) const
{
    /*vector<double> qq(15);
    for (auto i = 0; i < 15; ++i)
    {
        qq[i] = func(pts[i]);
    }
    const Point center = Point{ (pts[0].x + pts[1].x + pts[2].x) / 3, (pts[0].y + pts[1].y + pts[2].y) / 3 };
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0., 1.);
    double r1, r2;
    r1 = dist(mt);
    r2 = dist(mt);
    Point tp = center;
    //tp = Point((1 - sqrt(r1))*pts[0].x + sqrt(r1)*(1 - r2)*pts[1].x + sqrt(r1)*r2*pts[2].x, (1 - sqrt(r1))*pts[0].y + sqrt(r1)*(1 - r2)*pts[1].y + sqrt(r1)*r2*pts[2].y);
    //cout << GetShapeFunction(0, tp) << "\t" << GetShapeFunction(1, tp) << "\t" << GetShapeFunction(2, tp) << endl;
    //Point tp(Point{ pts[0].x + (pts[1].x - pts[0].x) / 2, pts[0].y + (pts[1].y - pts[0].y) / 2 });
    double expected = func(tp);
    Point dexpected(4. * tp.x * tp.x * tp.x, 4. * tp.y * tp.y * tp.y);
    double actual = 0;
    Point dactual;
    for (auto i = 0; i < m_s; ++i)
    {
        //cout << GetShapeFunction(i, tp) << endl;
        actual += qq[i] * GetShapeFunction(i, tp);
    }

    for (auto i = 0; i < m_s; ++i)
    {
        cout << GetGradShapeFunction(i, tp).x << "\t" << GetGradShapeFunction(i, tp).y << endl;
        dactual.x += qq[i] * GetGradShapeFunction(i, tp).x;
        dactual.y += qq[i] * GetGradShapeFunction(i, tp).y;
    }
    cout << node <<  "\t" << actual << "\t" << expected << endl;
    cout << dactual.x << "\t" << dexpected.x << endl;
    cout << dactual.y << "\t" << dexpected.y << endl;*/
    return func(pts[node]);
}
