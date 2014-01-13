#include <vector>

std::vector<std::vector<double> > convertRooTracker2DArray(double (*array)[4], int nRows, int nCols)
{
                    std::vector<std::vector<double> > result;
                    for(int i = 0; i < nRows; ++i)
                    {
                        std::vector<double> v;
                        for(int j = 0; j < nCols; ++j)
                        {
                            double d = array[i][j];
                            v.push_back(d);
                        }
                        result.push_back(v);
                    }
                    return result;
}

std::vector<std::vector<double> > getRooTrackerHepP4(int n, void* stdHepP4)
{
    int nRows = n;
    double (*array)[4] = (double (*)[4])stdHepP4;
    return convertRooTracker2DArray(array,nRows,4);
}
